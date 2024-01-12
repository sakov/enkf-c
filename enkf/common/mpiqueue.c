/******************************************************************************
 *
 * File:        mpiqueue.c        
 *
 * Created:     17/10/2023
 *
 * Author:      Pavel Sakov
 *
 * Description: Code for MPI job queue. Process with rank 0 manages the queue,
 *              the remaining CPUs do the jobs. See main() for example.
 *
 * Revisions:   
 *
 *****************************************************************************/

#if !defined(MPI)
typedef int make_iso_compilers_happy;
#else

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>
#include "mpiqueue.h"

#define MPIQUEUE_JOBSTATUS_TOASSIGN 0
#define MPIQUEUE_JOBSTATUS_ASSIGNED 1
#define MPIQUEUE_JOBSTATUS_DONE 2
#define MPIQUEUE_WORKERSTATUS_NA -1
#define MPIQUEUE_WORKERSTATUS_WAITING 0
#define MPIQUEUE_WORKERSTATUS_WORKING 1
#define MPIQUEUE_JOBTAG_OK 0
#define MPIQUEUE_JOBTAG_REJECT 1

struct mpiqueue {
    MPI_Comm communicator;
    int rank;
    int nprocesses;
    int* workerstatus;
    int njob;
    int* jobstatus;
    int mystatus;
};

static void quit_def(char* format, ...);
static mpiqueue_quit_fn quit = quit_def;

/**
 */
static void quit_def(char* format, ...)
{
    va_list args;

    fflush(stdout);

    fprintf(stderr, "\n\n  ERROR: mpiqueue: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n\n");
    fflush(NULL);
    exit(1);
}

/**
 */
mpiqueue* mpiqueue_create(MPI_Comm communicator, int njob)
{
    mpiqueue* queue = calloc(1, sizeof(mpiqueue));

    if (njob <= 0)
        quit("mpiqueue_create(): njob = %d; (should be > 0)", njob);
    queue->communicator = communicator;
    MPI_Comm_rank(communicator, &queue->rank);
    MPI_Comm_size(communicator, &queue->nprocesses);
    if (queue->nprocesses < 2)
        quit("mpiqueue_create(): nprocesses = %d; (should be > 1)", queue->nprocesses);
    queue->njob = njob;
    if (queue->rank == 0) {
        int jobid, p;

        queue->jobstatus = malloc(njob * sizeof(int));
        for (jobid = 0; jobid < njob; ++jobid)
            queue->jobstatus[jobid] = MPIQUEUE_JOBSTATUS_TOASSIGN;
        queue->mystatus = MPIQUEUE_WORKERSTATUS_NA;
        queue->workerstatus = malloc(queue->nprocesses * sizeof(int));
        queue->workerstatus[0] = MPIQUEUE_WORKERSTATUS_NA;
        for (p = 1; p < queue->nprocesses; ++p)
            queue->workerstatus[p] = MPIQUEUE_WORKERSTATUS_WAITING;
    } else {
        queue->mystatus = MPIQUEUE_WORKERSTATUS_WAITING;
    }

    /*
     * assign first batch of jobs
     */
    if (queue->rank == 0) {
        int jobid, p;

        for (p = 1, jobid = 0; p < queue->nprocesses && jobid < njob; ++p, ++jobid) {
            MPI_Send(&jobid, 1, MPI_INT, p, MPIQUEUE_JOBTAG_OK, queue->communicator);
            queue->jobstatus[jobid] = MPIQUEUE_JOBSTATUS_ASSIGNED;
            queue->workerstatus[p] = MPIQUEUE_WORKERSTATUS_WORKING;
        }
    }

    return queue;
}

/**
 */
void mpiqueue_destroy(mpiqueue* queue)
{
    if (queue->rank == 0) {
        free(queue->jobstatus);
        free(queue->workerstatus);
    }
    free(queue);
}

/**
 */
int mpiqueue_getrank(mpiqueue* queue)
{
    return queue->rank;
}

/**
 */
void mpiqueue_manage(mpiqueue* queue)
{
    if (queue->rank != 0)
        quit("mpiqueue_manage(): rank = %d (requires rank = 0)\n", queue->rank);

    while (1) {
        MPI_Status status;
        int jobid, j, p;

        /*
         * check whether all jobs are completed
         */
        for (jobid = 0; jobid < queue->njob; ++jobid)
            if (queue->jobstatus[jobid] != MPIQUEUE_JOBSTATUS_DONE)
                break;
        /*
         * if yes -- send finish signal to all workers and exit
         */
        if (jobid == queue->njob) {
            int finished = -1;

            /*
             * send completion signal
             */
            for (p = 1; p < queue->nprocesses; ++p)
                MPI_Send(&finished, 1, MPI_INT, p, 0, queue->communicator);
            return;
        }
        /*
         * if no -- wait for a message from workers
         */
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, queue->communicator, &status);
        MPI_Recv(&jobid, 1, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, queue->communicator, MPI_STATUS_IGNORE);
        if (jobid < 0 || jobid >= queue->njob)
            quit("jobid = %d (needs to be 0 <= jobid <= %d\n", jobid, queue->njob - 1);
        if (status.MPI_TAG == MPIQUEUE_JOBTAG_OK)
            queue->jobstatus[jobid] = MPIQUEUE_JOBSTATUS_DONE;
        else if (status.MPI_TAG == MPIQUEUE_JOBTAG_REJECT)
            queue->jobstatus[jobid] = MPIQUEUE_JOBSTATUS_TOASSIGN;
        queue->workerstatus[status.MPI_SOURCE] = MPIQUEUE_WORKERSTATUS_WAITING;

        jobid = (jobid + 1) % queue->njob;
        for (j = 0; j < queue->njob; ++j, jobid = (jobid + 1) % queue->njob)
            if (queue->jobstatus[jobid] == MPIQUEUE_JOBSTATUS_TOASSIGN)
                break;
        if (j == queue->njob)
            continue;

        /*
         * Because MPI_Probe() returned, there is at least one CPU waiting for
         * a new assignment. But we look for a spare CPU starting from the CPU
         * next to the one reported, to avoid locks when a CPU rejects a
         * certain job but gets it assigned over and over again.
         */
        /*
         * (skip the master -- CPU #0)
         */
        p = status.MPI_SOURCE % (queue->nprocesses - 1) + 1;
        for (j = 0; j < queue->nprocesses - 1; ++j, p = p % (queue->nprocesses - 1) + 1)
            if (queue->workerstatus[p] == MPIQUEUE_WORKERSTATUS_WAITING)
                break;
        {
            MPI_Request request;

            MPI_Isend(&jobid, 1, MPI_INT, p, 0, queue->communicator, &request);
            MPI_Request_free(&request);
        }
        queue->jobstatus[jobid] = MPIQUEUE_JOBSTATUS_ASSIGNED;
        queue->workerstatus[p] = MPIQUEUE_WORKERSTATUS_WORKING;
    }
}

/**
 */
int mpiqueue_getjob(mpiqueue* queue)
{
    int jobid;

    if (queue->rank == 0)
        quit("mpiqueue_getjobid(): called from process with rank = 0");
    if (queue->mystatus != MPIQUEUE_WORKERSTATUS_WAITING)
        quit("mpiqueue_getjobid(): worker #%d: requested new job without reporting completion of the previos one", queue->rank);
    MPI_Recv(&jobid, 1, MPI_INT, 0, MPI_ANY_TAG, queue->communicator, MPI_STATUS_IGNORE);
    queue->mystatus = MPIQUEUE_WORKERSTATUS_WORKING;

    return jobid;
}

/**
 */
void mpiqueue_reportjob(mpiqueue* queue, int jobid)
{
    if (queue->rank == 0)
        quit("mpiqueue_reportjobid(): called from master");
    if (queue->mystatus != MPIQUEUE_WORKERSTATUS_WORKING)
        quit("mpiqueue_reportjobid(): worker #%d: reporting completion of job #%d without being assigned", queue->rank, jobid);
    MPI_Send(&jobid, 1, MPI_INT, 0, 0, queue->communicator);
    queue->mystatus = MPIQUEUE_WORKERSTATUS_WAITING;
}

/**
 */
void mpiqueue_rejectjob(mpiqueue* queue, int jobid)
{
    if (queue->rank == 0)
        quit("mpiqueue_reportjobid(): called from master");
    MPI_Send(&jobid, 1, MPI_INT, 0, MPIQUEUE_JOBTAG_REJECT, queue->communicator);
    queue->mystatus = MPIQUEUE_WORKERSTATUS_WAITING;
}

/**
 */
void mpiqueue_setquitfn(mpiqueue_quit_fn quit_fn)
{
    quit = quit_fn;
}

#if defined MPIQUEUE_TEST
/**
 */
#include <limits.h>
#include <unistd.h>

/**
 */
static int str2int(char* token, int* value)
{
    long int tmp;
    char* end = NULL;

    if (token == NULL) {
        *value = INT_MAX;
        return 0;
    }

    tmp = strtol(token, &end, 10);

    if (end == token || tmp > INT_MAX || tmp < INT_MIN) {
        *value = INT_MAX;
        return 0;
    }

    *value = (int) tmp;
    return 1;
}

/**
 */
static void dojob(int rank, int jobid, int taskid)
{
    size_t answer;
    int i;

    if (jobid <= 20) {
        answer = 1;
        for (i = 2; i <= jobid; ++i)
            answer *= i;
        printf("  task %d: I am #%d: %d! = %zu\n", taskid, rank, jobid, answer);
    } else
        printf("  task %d: I am #%d: %d! = N/A (overflow)\n", taskid, rank, jobid);

    fflush(stdout);
}

/**
 */
int main(int argc, char* argv[])
{
    int N;
    mpiqueue* queue = NULL;
    int rank;

    if (argc != 2) {
        printf("  Usage: mpirun -np <NCPU> mpiqueue_test <N>\n");
        exit(argc != 1);
    }
    if (!str2int(argv[1], &N))
        quit("could not convert \"%s\" to int", argv[1]);

    MPI_Init(&argc, &argv);

    /*
     * jobs: calculate n! for n = 0 ... N
     */
    queue = mpiqueue_create(MPI_COMM_WORLD, N + 1);
    rank = mpiqueue_getrank(queue);
    if (rank == 0) {
        printf("  task 1:\n");
        printf("    calculate n! for n = 0...%d\n", N);
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    /*
     * main cycle
     */
    if (rank == 0)
        mpiqueue_manage(queue);
    else {
        int jobid;

        while ((jobid = mpiqueue_getjob(queue)) >= 0) {
            dojob(rank, jobid, 1);
            mpiqueue_reportjob(queue, jobid);
        }
    }
    mpiqueue_destroy(queue);

    MPI_Barrier(MPI_COMM_WORLD);
    usleep(100);

    /*
     * jobs: calculate n! for n = 0...N
     * but only allow even CPUs calculate factorial for even n,
     * odd CPUs calculate factorial for odd n
     */
    queue = mpiqueue_create(MPI_COMM_WORLD, N + 1);
    if (queue->nprocesses < 3)
        goto finish;
    rank = mpiqueue_getrank(queue);
    if (rank == 0) {
        printf("\n  task 2:\n");
        printf("    calculate n! for n = 0...%d\n", N);
        printf("      allow CPUs with odd rank work on odd n only\n");
        printf("      allow CPUs with even rank work on even n only\n");
        fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        mpiqueue_manage(queue);
    else {
        int jobid;

        while ((jobid = mpiqueue_getjob(queue)) >= 0) {
            if (jobid % 2 != rank % 2) {
                mpiqueue_rejectjob(queue, jobid);
                continue;
            }
            dojob(rank, jobid, 2);
            mpiqueue_reportjob(queue, jobid);
        }
    }
  finish:
    mpiqueue_destroy(queue);

    MPI_Finalize();

    return 0;
}
#endif                          /* MPIQUEUE_TEST */
#endif                          /* MPI */
