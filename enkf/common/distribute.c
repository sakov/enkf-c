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

void distribute_iterations(int i1, int i2, int nproc, int rank, char prefix[])
{
    int n, npp, i, j;

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
    if (prefix != NULL)
        enkf_printf("%sdistributing iterations:\n", prefix);
#if defined(MPI)
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    n = i2 - i1 + 1;
    npp = n / nproc;
    j = -1;
    if (n % nproc == 0) {
        my_number_of_iterations = n / nproc;
        for (i = 0; i < nproc; ++i)
            number_of_iterations[i] = my_number_of_iterations;
        if (prefix != NULL)
            enkf_printf("%s  all processes get %d iterations\n", prefix, my_number_of_iterations);
        j = nproc;
    } else {
        for (i = 1; i < nproc; ++i) {
            if (i * (npp + 1) + (nproc - i) * npp == n) {
                j = i;
                break;
            }
        }

        if (rank < j)
            my_number_of_iterations = npp + 1;
        else
            my_number_of_iterations = npp;

        for (i = 0; i < j; ++i)
            number_of_iterations[i] = npp + 1;
        for (i = j; i < nproc; ++i)
            number_of_iterations[i] = npp;
        assert(j * (npp + 1) + (nproc - j) * npp == n);
        if (prefix != NULL)
            enkf_printf("%s  processes get %d or %d iterations\n", prefix, number_of_iterations[0], number_of_iterations[nproc - 1]);
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
    number_of_iterations = NULL;
    free(first_iteration);
    first_iteration = NULL;
    free(last_iteration);
    last_iteration = NULL;
    my_number_of_iterations = -1;
    my_first_iteration = -1;
    my_last_iteration = -1;
}
