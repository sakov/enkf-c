/******************************************************************************
 *
 * File:        mpiqueue.h        
 *
 * Created:     17/10/2023
 *
 * Author:      Pavel Sakov
 *
 * Description: Code for MPI job queue. Process with rank 0 manages the queue,
 *              the remaining CPUs do the jobs. See main() in mpiqueue.c for
 *              example.
 *
 * Revisions:   
 *
 *****************************************************************************/

#if !defined(_MPIQUEUE_H)
#define _MPIQUEUE_H

struct mpiqueue;
typedef struct mpiqueue mpiqueue;

#if defined(__GNUC__)
typedef void (*mpiqueue_quit_fn) (char* format, ...) __attribute__ ((format(printf, 1, 2)));
#else
typedef void (*mpiqueue_quit_fn) (char* format, ...);
#endif
void mpiqueue_setquitfn(mpiqueue_quit_fn quit_fn);

mpiqueue* mpiqueue_create(MPI_Comm communicator, int njob);
void mpiqueue_destroy(mpiqueue* queue);
int mpiqueue_getrank(mpiqueue* queue);
void mpiqueue_manage(mpiqueue* queue);
int mpiqueue_getjob(mpiqueue* queue);
void mpiqueue_reportjob(mpiqueue* queue, int jobid);
void mpiqueue_rejectjob(mpiqueue* queue, int jobid);

#endif                          /* _MPIQUEUEH */
