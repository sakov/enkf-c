/******************************************************************************
 *
 * File:        distribute.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: A simple procedure for distributing iterations between
 *              processes. It assumes that
 *              (1) there is only one job handled by this module at a time;
 *              (2) the process ID and the total number of available slots are
 *                  stored in global variables `rank' and `nprocesses'.
 *
 *              The iterations in the interval [i1, i2] are distributed between
 *              `nslot' units (CPUs) (nslot <= nprocesses) with IDs
 *               0 <= rank < nslot.
 *
 *              The results are stored in 6 global variables, with the following
 *              relations between them:
 *                my_number_of_iterations = my_last_iteration 
 *                                                      - my_first_iteration + 1
 *                my_number_of_iterations = number_of_iterations[rank]
 *                my_first_iteration = first_iteration[rank]
 *                my_last_iteration = last_iteration[rank]
 *
 * Revisions:   18/04/2018 PS: `nproc' was supposed to be an alias for
 *              `nprocesses'; now it can be arbitrary number such that
 *              0 < nproc < nprocesses.
 *              21/11/2024 PS: use `nslot' instead of `nproc'.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <assert.h>
#if defined(MPI)
#include <mpi.h>
#endif
#include "utils.h"

extern int nprocesses;          /* the total number of available slots */
extern int rank;                /* own rank */

int my_number_of_iterations = -1;
int my_first_iteration = -1;
int my_last_iteration = -1;
int* number_of_iterations = NULL;
int* first_iteration = NULL;
int* last_iteration = NULL;

/** Distributes indices in the interval [i1, i2] between slots in the interval
 ** [0, nslot - 1]. The results are written to 6 global variables.
 * @param i1 Start of the interval
 * @param i2 End of the interval
 * @param nslot Number of slots (CPUs) to be be used
 * @param prefix Prefix for log printing; NULL to print no log.
 * There also need to be two global variables:
 *   nprocesses -- the total number of available slots;
 *   rank -- the ID of own slot.
 */
void distribute_iterations(int i1, int i2, int nslot, char prefix[])
{
    int niter = i2 - i1 + 1;
    int ii, i, jj;

    assert(i2 >= i1);
    assert(nslot > 0 && nslot <= nprocesses);

#if defined(MPI)
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if (number_of_iterations == NULL) {
        number_of_iterations = malloc(nprocesses * sizeof(int));
        first_iteration = malloc(nprocesses * sizeof(int));
        last_iteration = malloc(nprocesses * sizeof(int));
    }

    if (prefix != NULL)
        enkf_printf("%sdistributing %d iterations:\n", prefix, niter);
#if defined(MPI)
    if (prefix != NULL)
        fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    ii = niter / nslot;
    jj = niter % nslot;
    for (i = 0; i < jj; ++i)
        number_of_iterations[i] = ii + 1;
    for (i = jj; i < nslot; ++i)
        number_of_iterations[i] = ii;
    for (i = nslot; i < nprocesses; ++i)
        number_of_iterations[i] = 0;

    first_iteration[0] = i1;
    last_iteration[0] = i1 + number_of_iterations[0] - 1;
    for (i = 1; i < nprocesses; ++i) {
        first_iteration[i] = last_iteration[i - 1] + 1;
        last_iteration[i] = first_iteration[i] + number_of_iterations[i] - 1;
    }

    my_first_iteration = first_iteration[rank];
    my_last_iteration = last_iteration[rank];
    my_number_of_iterations = number_of_iterations[rank];
    
    if (prefix != NULL) {
        if (nslot == 1)
            enkf_printf("%s  1 process gets %d iterations\n", prefix, ii);
        else if (jj == 0)
            enkf_printf("%s  %d processes get %d iterations each\n", prefix, nslot, ii);
        else if (ii == 0)
            enkf_printf("%s  %d processes get 1 iteration, %d processes get 0 iterations\n", prefix, niter, nslot - niter);
        else
            enkf_printf("%s  %d processes get %d or %d iterations\n", prefix, nslot, ii + 1, ii);
        if (nprocesses > nslot)
            enkf_printf("%s  %d processes not involved\n", prefix, nprocesses - nslot);
        enkf_flush();
    }
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

/**
 */
void distribute_free(void)
{
    if (number_of_iterations == NULL)
        return;
    free(number_of_iterations);
    number_of_iterations = NULL;
    free(first_iteration);
    first_iteration = NULL;
    free(last_iteration);
    last_iteration = NULL;
    my_number_of_iterations = -1;
    my_first_iteration = -1;
    my_last_iteration = -1;
}
