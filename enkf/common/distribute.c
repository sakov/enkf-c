/******************************************************************************
 *
 * File:        distribute.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
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

void distribute_iterations(int i1, int i2, int nproc, int rank)
{
    double nproc_dbl, npp;
    int n, i, j;

#if defined(MPI)
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    assert(i2 >= i1);

    if (number_of_iterations == NULL) {
        number_of_iterations = malloc(nproc * sizeof(int));
        first_iteration = malloc(nproc * sizeof(int));
        last_iteration = malloc(nproc * sizeof(int));
    }
    if (rank == 0)
        enkf_printf("    distributing iterations:\n");
#if defined(MPI)
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    nproc_dbl = (double) nproc;
    n = i2 - i1 + 1;
    npp = (double) n / nproc_dbl;

    j = -1;
    if (floor(npp) == npp) {
        my_number_of_iterations = n / nproc;
        for (i = 0; i < nproc; ++i)
            number_of_iterations[i] = my_number_of_iterations;
        enkf_printf("      all processes get %d iterations\n", my_number_of_iterations);
        j = nproc;
    } else {
        for (i = 1; i < nproc; ++i) {
            if (i * ((int) npp + 1) + (nproc - i) * (int) npp == n) {
                j = i;
                break;
            }
        }

        if (rank < j)
            my_number_of_iterations = (int) npp + 1;
        else
            my_number_of_iterations = (int) npp;

        for (i = 0; i < j; ++i)
            number_of_iterations[i] = (int) npp + 1;
        for (i = j; i < nproc; ++i)
            number_of_iterations[i] = (int) npp;
        assert(j * ((int) npp + 1) + (nproc - j) * (int) npp == n);
        if (n < nproc)
            enkf_quit("distribute_iterations(): number of elements less than number of processors");
    }
#if defined(MPI)
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    first_iteration[0] = i1;
    last_iteration[0] = i1 + number_of_iterations[0] - 1;
    for (i = 1; i < nproc; ++i) {
        first_iteration[i] = last_iteration[i - 1] + 1;
        last_iteration[i] = first_iteration[i] + number_of_iterations[i] - 1;
    }

    my_first_iteration = first_iteration[rank];
    my_last_iteration = last_iteration[rank];

    enkf_printf("      process %d: %d - %d\n", rank, my_first_iteration, my_last_iteration);

#if defined(MPI)
    fflush(stdout);
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
    free(first_iteration);
    free(last_iteration);
}
