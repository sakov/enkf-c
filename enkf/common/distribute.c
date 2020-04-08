/******************************************************************************
 *
 * File:        distribute.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Distributes indices in the interval [i1, i2] between `nproc'
 *              processes. Process IDs are assumed to be in the interval
 *              [0, nproc-1]. The calling process has ID `rank'. The results
 *              are stored in 6 global variables, with the following relations
 *              between them:
 *                my_number_of_iterations = my_last_iteration 
 *                                                      - my_first_iteration + 1
 *                my_number_of_iterations = number_of_iterations[rank]
 *                my_first_iteration = first_iteratin[rank]
 *                my_last_iteration = last_iteration[rank]
 *
 * Revisions:   18/04/2018 PS: `nproc' was supposed to be an alias for
 *              `nprocesses'; now it can be arbitrary number such that
 *              0 < nproc < nprocesses.
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
#include "definitions.h"
#include "utils.h"

int my_number_of_iterations = -1;
int my_first_iteration = -1;
int my_last_iteration = -1;
int* number_of_iterations = NULL;
int* first_iteration = NULL;
int* last_iteration = NULL;

/** Distributes indices in the interval [i1, i2] between `nproc' processes.
 * @param i1 Start of the interval
 * @param i2 End of the interval
 * @param nproc Number of processes (CPUs) to be be used
 * @param rank ID of the "native" process
 * @param prefix Prefix for log printing; NULL to print no log.
 * Note that `nprocesses' and `rank' are supposed to be external (global) 
 * variables.
 */
void distribute_iterations(int i1, int i2, int nproc, int rank, char prefix[])
{
    int n, npp, i;

#if defined(MPI)
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    assert(i2 >= i1);

    if (number_of_iterations == NULL) {
        number_of_iterations = malloc(nprocesses * sizeof(int));
        first_iteration = malloc(nprocesses * sizeof(int));
        last_iteration = malloc(nprocesses * sizeof(int));
    }
    if (prefix != NULL)
        enkf_printf("%sdistributing %d iterations:\n", prefix, i2 - i1 + 1);
#if defined(MPI)
    if (prefix != NULL)
        fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    assert(nproc > 0 && nproc <= nprocesses);

    n = i2 - i1 + 1;
    npp = n / nproc;
    if (n % nproc == 0) {
        for (i = 0; i < nproc; ++i)
            number_of_iterations[i] = npp;
        if (prefix != NULL)
            enkf_printf("%s  %d processes get %d iterations each\n", prefix, nproc, npp);
    } else {
        int j = n - nproc * npp;

        for (i = 0; i < j; ++i)
            number_of_iterations[i] = npp + 1;
        for (i = j; i < nproc; ++i)
            number_of_iterations[i] = npp;
        if (prefix != NULL) {
            if (npp == 0) {
                enkf_printf("%s  %d processes get 1 iteration, %d processes get 0 iterations\n", prefix, n, nproc - n);
            } else
                enkf_printf("%s  %d processes get %d or %d iterations\n", prefix, nproc, npp + 1, npp);
        }
    }
    for (i = nproc; i < nprocesses; ++i)
        number_of_iterations[i] = 0;
    if (nprocesses > nproc && prefix != NULL)
        enkf_printf("%s  %d processes get 0 iterations\n", prefix, nprocesses - nproc);
#if defined(MPI)
    if (prefix != NULL)
        enkf_flush();
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    first_iteration[0] = i1;
    last_iteration[0] = i1 + number_of_iterations[0] - 1;
    for (i = 1; i < nproc; ++i) {
        first_iteration[i] = last_iteration[i - 1] + 1;
        last_iteration[i] = first_iteration[i] + number_of_iterations[i] - 1;
    }
    for (i = nproc; i < nprocesses; ++i) {
        first_iteration[i] = last_iteration[i - 1] + 1;
        last_iteration[i] = first_iteration[i] + number_of_iterations[i] - 1;
    }

    my_first_iteration = first_iteration[rank];
    my_last_iteration = last_iteration[rank];
    my_number_of_iterations = number_of_iterations[rank];
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
