/******************************************************************************
 *
 * File:        mpiqueue.c        
 *
 * Created:     17/10/2023
 *
 * Author:      Pavel Sakov
 *
 * Description: Code for MPI job queue. Process with rank 0 manages the queue,
 *              other processes do the jobs. See main() for example.
 *
 * Revisions:   
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>
#include "mpiqueue.h"

#define MPIQUEUE_JOBSTATUS_TOASSIGN 0
#define MPIQUEUE_JOBSTATUS_ASSIGNED 1
#define MPIQUEUE_JOBSTATUS_DONE 2

struct mpiqueue {
    MPI_Comm communicator;
    int rank;
    int nprocesses;
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
    mpiqueue* queue = malloc(sizeof(mpiqueue));

    if (njob <= 0)
        quit("mpiqueue_create(): njob = %d; (should be > 0)", njob);
    queue->communicator = communicator;
    MPI_Comm_rank(communicator, &queue->rank);
    MPI_Comm_size(communicator, &queue->nprocesses);
    if (queue->nprocesses < 2)
        quit("nprocesses = %d; (should be > 1)", queue->nprocesses);
    queue->njob = njob;
    if (queue->rank == 0) {
        int jobid;

        queue->jobstatus = malloc(njob * sizeof(int));
        for (jobid = 0; jobid < njob; ++jobid)
            queue->jobstatus[jobid] = MPIQUEUE_JOBSTATUS_TOASSIGN;
    } else
        queue->jobstatus = NULL;

    /*
     * assign first batch of jobs
     */
    if (queue->rank == 0) {
        int jobid, p;

        for (p = 1, jobid = 0; p < queue->nprocesses && jobid < njob; ++p, ++jobid) {
            MPI_Send(&jobid, 1, MPI_INT, p, 0, queue->communicator);
            queue->jobstatus[jobid] = MPIQUEUE_JOBSTATUS_ASSIGNED;
        }
    }

    return queue;
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
        int jobid;

        for (jobid = 0; jobid < queue->njob; ++jobid)
            if (queue->jobstatus[jobid] != MPIQUEUE_JOBSTATUS_DONE)
                break;
        if (jobid == queue->njob) {
            int finished = -1;
            int p;

            for (p = 1; p < queue->nprocesses; ++p)
                MPI_Send(&finished, 1, MPI_INT, p, 0, queue->communicator);
            return;
        }

        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, queue->communicator, &status);
        MPI_Recv(&jobid, 1, MPI_INT, status.MPI_SOURCE, MPI_ANY_TAG, queue->communicator, MPI_STATUS_IGNORE);
        if (jobid < 0 || jobid >= queue->njob)
            quit("jobid = %d (needs to be 0 <= jobid <= %d\n", jobid, queue->njob - 1);
        queue->jobstatus[jobid] = MPIQUEUE_JOBSTATUS_DONE;

        for (jobid = 0; jobid < queue->njob; ++jobid)
            if (queue->jobstatus[jobid] == MPIQUEUE_JOBSTATUS_TOASSIGN)
                break;
        if (jobid == queue->njob)
            continue;
        MPI_Send(&jobid, 1, MPI_INT, status.MPI_SOURCE, 0, queue->communicator);
        queue->jobstatus[jobid] = MPIQUEUE_JOBSTATUS_ASSIGNED;
    }
}

/**
 */
int mpiqueue_getjobid(mpiqueue* queue)
{
    int jobid;

    if (queue->rank == 0)
        quit("mpiqueue_getjobid(): called from process with rank = 0");
    MPI_Recv(&jobid, 1, MPI_INT, 0, MPI_ANY_TAG, queue->communicator, MPI_STATUS_IGNORE);

    return jobid;
}

/**
 */
void mpiqueue_reportjobid(mpiqueue* queue, int jobid)
{
    if (queue->rank == 0)
        quit("mpiqueue_reportjobid(): called from process with rank = 0");
    MPI_Send(&jobid, 1, MPI_INT, 0, 0, queue->communicator);
}

/**
 */
void mpiqueue_destroy(mpiqueue* queue)
{
    if (queue->rank == 0)
        free(queue->jobstatus);
    free(queue);
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
int main(int argc, char* argv[])
{
    int N;
    mpiqueue* queue = NULL;

    if (argc != 2) {
        printf("  Usage: mpirun -np <NCPU> mpiqueue_test <N>\n");
        exit(argc != 1);
    }
    if (!str2int(argv[1], &N))
        quit("could not convert \"%s\" to int", argv[1]);

    MPI_Init(&argc, &argv);
    queue = mpiqueue_create(MPI_COMM_WORLD, N);

    if (mpiqueue_getrank(queue) == 0)
        mpiqueue_manage(queue);
    else {
        int jobid;

        while ((jobid = mpiqueue_getjobid(queue)) >= 0) {
            size_t answer;
            int i;

            answer = 1;
            for (i = 2; i <= jobid; ++i)
                answer *= i;
            printf("  I am #%d: %d! = %zu\n", queue->rank, jobid, answer);
            fflush(stdout);
            mpiqueue_reportjobid(queue, jobid);
        }
    }

    mpiqueue_destroy(queue);
    MPI_Finalize();

    return 0;
}
#endif                          /* MPIQUEUE_TEST */
