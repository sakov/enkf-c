/******************************************************************************
 *
 * File:        kdtree.c        
 *
 * Created:     12/2012
 *
 * Author:      John Tsiombikas <nuclear@siggraph.org>
 *
 * Description:
 *
 * Revisions:   April 2012; March 2016 -- modified by Pavel Sakov
 *
 *****************************************************************************/

/*
This file is part of ``kdtree'', a library for working with kd-trees.
Copyright (C) 2007-2009 John Tsiombikas <nuclear@siggraph.org>

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
OF SUCH DAMAGE.
*/
/* single nearest neighbor search written by Tamas Nepusz <tamas@cs.rhul.ac.uk> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kdtree.h"

#define NALLOCSTART 1024

struct resnode;
typedef struct resnode resnode;

struct kdnode {
    int id;
    int dir;
    int left;
    int right;
};

struct kdtree {
    int ndim;
    int nnodes;
    kdnode* nodes;
    double* positions;
    int nallocated;
};

struct resnode {
    int id;
    double dist_sq;
    resnode* next;
};

struct kdset {
    struct kdtree* tree;
    resnode* first, *now;
    int size;
};

/**
 */
kdtree* kd_create(int ndim)
{
    kdtree* tree;

    tree = malloc(sizeof(kdtree));
    tree->ndim = ndim;
    tree->nnodes = 0;
    tree->nodes = NULL;
    tree->positions = NULL;
    tree->nallocated = 0;

    return tree;
}

/**
 */
void kd_destroy(kdtree* tree)
{
    free(tree->nodes);
    free(tree->positions);
    free(tree);
}

/**
 */
static int _kd_insert(kdtree* tree, int* id, const double* pos, int dir)
{
    int ndim = tree->ndim;
    int new_dir;
    kdnode* node;

    if (*id == -1) {
        node = &tree->nodes[tree->nnodes];
        memcpy(&tree->positions[tree->nnodes * ndim], pos, ndim * sizeof(double));
        node->id = tree->nnodes;
        node->dir = dir;
        node->left = -1;
        node->right = -1;
        *id = node->id;
        return 0;
    }
    node = &tree->nodes[*id];
    new_dir = (node->dir + 1) % ndim;

    return _kd_insert(tree, (pos[node->dir] < tree->positions[node->id * ndim + node->dir]) ? &node->left : &node->right, pos, new_dir);
}

/**
 */
int kd_insert(kdtree* tree, const double* pos)
{
    int status;

    if (tree->nallocated == tree->nnodes) {
        if (tree->nallocated != 0)
            tree->nallocated *= 2;
        else
            tree->nallocated = NALLOCSTART;
        tree->nodes = realloc(tree->nodes, tree->nallocated * sizeof(kdnode));
        tree->positions = realloc(tree->positions, tree->nallocated * tree->ndim * sizeof(double));
        if (tree->nallocated == NALLOCSTART)
            tree->nodes[0].id = -1;
    }
    status = _kd_insert(tree, &tree->nodes[0].id, pos, 0);
    if (status == 0)
        tree->nnodes++;

    return status;
}

/** insert the item. if dist_sq is >= 0, then do an ordered insert
 */
static int _rlist_insert(resnode* list, kdnode* item, double dist_sq)
{
    resnode* rnode = malloc(sizeof(resnode));

    rnode->id = item->id;
    rnode->dist_sq = dist_sq;

    if (dist_sq >= 0.0)
        while (list->next && list->next->dist_sq < dist_sq)
            list = list->next;

    rnode->next = list->next;
    list->next = rnode;

    return 0;
}

/**
 */
static int _kd_nearest_range(kdtree* tree, int id, const double* pos, double range, resnode* list, int ordered)
{
    int ndim = tree->ndim;
    kdnode* node;
    double* nodepos;
    double dist_sq, dx;
    int i, ret, added_res;

    if (id < 0)
        return 0;

    node = &tree->nodes[id];
    nodepos = &tree->positions[node->id * ndim];

    dist_sq = 0;
    for (i = 0; i < ndim; i++) {
        double tmp = nodepos[i] - pos[i];

        dist_sq += tmp * tmp;
    }

    added_res = 0;
    if (dist_sq <= range * range) {
        if (_rlist_insert(list, node, ordered ? dist_sq : -1.0) == -1)
            return -1;
        added_res = 1;
    }

    dx = pos[node->dir] - nodepos[node->dir];
    ret = _kd_nearest_range(tree, dx <= 0.0 ? node->left : node->right, pos, range, list, ordered);
    if (ret >= 0 && fabs(dx) < range) {
        added_res += ret;
        ret = _kd_nearest_range(tree, dx <= 0.0 ? node->right : node->left, pos, range, list, ordered);
    }
    if (ret == -1)
        return -1;
    added_res += ret;

    return added_res;
}

/**
 */
kdset* kd_nearest_range(kdtree* tree, const double* pos, double range, int ordered)
{
    int ret;
    kdset* rset;

    rset = malloc(sizeof(kdset));
    rset->first = malloc(sizeof(resnode));
    rset->first->next = NULL;
    rset->tree = tree;

    ret = _kd_nearest_range(tree, 0, pos, range, rset->first, ordered);
    if (ret == -1) {
        kd_res_free(rset);
        return 0;
    }
    rset->size = ret;
    rset->now = rset->first->next;    /* rewind */

    return rset;
}

/**
 */
static void _clear_results(kdset* rset)
{
    resnode* tmp, *node = rset->first->next;

    while (node != NULL) {
        tmp = node;
        node = node->next;
        free(tmp);
    }
    rset->first->next = NULL;
}

/**
 */
void kd_res_free(kdset* set)
{
    if (set == NULL)
        return;
    _clear_results(set);
    free(set->first);
    free(set);
}

/**
 */
int kd_res_next(kdset* set)
{
    if (set == NULL || set->now == NULL)
        return 0;
    set->now = set->now->next;

    return set->now != NULL;
}

/**
 */
int kd_res_getid(kdset* set)
{
    return (set->now != NULL) ? set->now->id : -1;
}

/**
 */
int kd_getid(kdnode* node)
{
    return node->id;
}

/**
 */
double* kd_getpos(kdtree* tree, int id)
{
    return &tree->positions[id * tree->ndim];
}

#if defined(STANDALONE)

#include <stdarg.h>
#include <errno.h>
#include "nan.h"

#define BUFSIZE 1024
#define NDIMMAX 256
#define NDATAMAX 256

#define INVRAD (3.14159265358979323846 / 180.0)
#define REARTH 6371.0

static int isll = 0;

/**
 */
static void quit(char* format, ...)
{
    va_list args;

    fflush(stdout);

    fprintf(stderr, "  error: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n\n");
    exit(1);
}

/**
 */
static int str2double(char* token, double* value)
{
    char* end = NULL;

    if (token == NULL) {
        *value = NaN;
        return 0;
    }

    *value = strtod(token, &end);

    if (end == token) {
        *value = NaN;
        return 0;
    }

    return 1;
}

/**
 */
static void* kd_readfile(char* fname, int ndim)
{
    FILE* f = NULL;
    kdtree* tree = NULL;
    char buf[BUFSIZE];
    char seps[] = " ,;\t";
    char* token;
    double point[NDIMMAX];
    int i;

    if (ndim < 1 || ndim > NDIMMAX)
        quit("no. of dimensions = %d; expected 1 <= ndim <= %d", ndim, NDIMMAX);
    if (fname == NULL)
        f = stdin;
    else {
        if (strcmp(fname, "stdin") == 0 || strcmp(fname, "-") == 0)
            f = stdin;
        else {
            f = fopen(fname, "r");
            if (f == NULL)
                quit("%s: %s", fname, strerror(errno));
        }
    }

    tree = kd_create((isll) ? 3 : ndim);
    while (fgets(buf, BUFSIZE, f) != NULL) {
        if (buf[0] == '#')
            continue;

        for (i = 0; i < ndim; ++i) {
            if ((token = strtok((i == 0) ? buf : NULL, seps)) == NULL)
                break;
            if (!str2double(token, &point[i]))
                break;
        }
        if (i < ndim)
            continue;

        if (isll) {
            double lon = point[0] * INVRAD;
            double lat = point[1] * INVRAD;
            double coslat = cos(lat);

            point[0] = REARTH * sin(lon) * coslat;
            point[1] = REARTH * cos(lon) * coslat;
            point[2] = REARTH * sin(lat);
        }

        kd_insert(tree, point);
    }

    if (f != stdin)
        if (fclose(f) != 0)
            quit("%s: %s", fname, strerror(errno));

    return tree;
}

/**
 */
static void description(void)
{
    printf("  kd: reads coordinates from a text file; searches for points\n");
    printf("      within specified radius from a specified location; prints results\n");
    printf("      (point coordinates and id) to the standard output\n");
}

/**
 */
static void usage(void)
{
    printf("  Usage: kd -i <file> -n <ndim> -p <coords> -r <range> [-s] [-g]\n");
    printf("         kd -h\n");
    printf("  Options:\n");
    printf("    -i <file>   -- a text file, with each row containing <ndim> point\n");
    printf("                   coordinates\n");
    printf("    -g          -- assume coordinates to be lon/lat in degrees and the distance is in kilometres\n");
    printf("    -h          -- describe this program\n");
    printf("    -n <ndim>   -- number of dimensions\n");
    printf("    -p <coords> -- <ndim> coordinates of the point to search around\n");
    printf("    -r <range>  -- radius to search within\n");
    printf("    -s          -- sort results\n");
    printf("                   <range> in km\n");

    exit(0);
}

/**
 */
static void parse_commandline(int argc, char* argv[], char** fname, int* ndim, double* pos, double* range, int* dosort)
{
    int i, ii;

    if (argc < 2)
        usage();

    i = 1;
    while (i < argc) {
        if (argv[i][0] != '-') {
            printf("\n%s: error for argument # %d\n\n", argv[0], i);
            usage();
        }
        switch (argv[i][1]) {
        case 'i':
            i++;
            if (i >= argc)
                quit("no <fname> value found after -i");
            *fname = argv[i];
            i++;
            break;
        case 'n':
            i++;
            if (i >= argc)
                quit("no <ndim> value found after -n");
            *ndim = atoi(argv[i]);
            i++;
            break;
        case 'p':
            i++;
            if (i >= argc)
                quit("no <coords> values found after -p");
            for (ii = 0; argv[i][0] != '-'; ++ii, ++i)
                if (!str2double(argv[i], &pos[ii]))
                    usage();
            break;
        case 'r':
            i++;
            if (i >= argc)
                quit("no <range> value found after -d");
            if (!str2double(argv[i], range))
                usage();
            i++;
            break;
        case 's':
            i++;
            *dosort = 1;
            break;
        case 'g':
            i++;
            isll = 1;
            if (*ndim <= 0)
                *ndim = 2;
            break;
        case 'h':
            description();
            exit(0);
        default:
            usage();
        }
    }

    if (*ndim <= 0)
        quit("\"ndim\" must be defined");

    if (!isfinite(*range) || *range <= 0.0)
        quit("range invalid or not defined");
}

/**
 */
int main(int argc, char* argv[])
{
    char* fname = NULL;
    int ndim = -1;
    double pos[NDIMMAX];
    double range = NaN;
    int dosort = 0;
    kdtree* tree = NULL;
    kdset* set = NULL;
    int id, i;

    parse_commandline(argc, argv, &fname, &ndim, pos, &range, &dosort);

    tree = kd_readfile(fname, ndim);

    if (isll) {
        double lon = pos[0] * INVRAD;
        double lat = pos[1] * INVRAD;

        pos[0] = REARTH * sin(lon) * cos(lat);
        pos[1] = REARTH * cos(lon) * cos(lat);
        pos[2] = REARTH * sin(lat);
    }

    set = kd_nearest_range(tree, pos, range, dosort);

    for (; (id = kd_res_getid(set)) >= 0; kd_res_next(set)) {
        double* respos;
        double dist;

        respos = kd_getpos(tree, id);

        if (isll) {
            double lat = asin(respos[2] / REARTH);
            double lon = asin(respos[0] / REARTH / cos(lat));

            dist = sqrt((respos[0] - pos[0]) * (respos[0] - pos[0]) + (respos[1] - pos[1]) * (respos[1] - pos[1]) + (respos[2] - pos[2]) * (respos[2] - pos[2]));

            respos[0] = lon / INVRAD;
            respos[1] = lat / INVRAD;
        } else {
            double dist_sq = 0.0;

            for (i = 0; i < ndim; ++i)
                dist_sq += (respos[i] - pos[i]) * (respos[i] - pos[i]);
            dist = sqrt(dist_sq);
        }

        for (i = 0; i < ndim; ++i)
            printf("%10f ", respos[i]);
        printf("%4d ", id);
        printf("%g\n", dist);
    }

    kd_res_free(set);
    kd_destroy(tree);

    return 0;
}
#endif                          /* STANDALONE */
