/******************************************************************************
 *
 * File:        mpiqueue.h        
 *
 * Created:     17/10/2023
 *
 * Author:      Pavel Sakov
 *
 * Description: Code for MPI job queue. Process with rank 0 manages the queue,
 *              other processes do the jobs. See main() in mpiqueue.c for
 *              example.
 *
 * Revisions:   
 *
 *****************************************************************************/

#if !defined(_MPIQUEUE_H)
#define _MPIQUEUE_H

struct mpiqueue;
typedef struct mpiqueue mpiqueue;

typedef void (*mpiqueue_quit_fn) (char* format, ...);
void mpiqueue_setquitfn(mpiqueue_quit_fn quit_fn);

mpiqueue* mpiqueue_create(MPI_Comm communicator, int njob);
int mpiqueue_getrank(mpiqueue* queue);
void mpiqueue_manage(mpiqueue* queue);
int mpiqueue_getjobid(mpiqueue* queue);
void mpiqueue_reportjobid(mpiqueue* queue, int jobid);
void mpiqueue_destroy(mpiqueue* queue);

#endif                          /* _MPIQUEUEH */
