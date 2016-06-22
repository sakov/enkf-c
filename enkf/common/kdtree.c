/******************************************************************************
 *
 * File:        kdtree.c        
 *
 * Created:     23/03/2016
 *
 * Author:      Pavel Sakov
 *              Derived from the code by John Tsiombikas (see the tail of the
 *              file)
 *
 * Description: KD-tree code
 *
 * Revisions:   
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include "kdtree.h"

#define NALLOCSTART 1024
#define SEED 5555

struct resnode;
typedef struct resnode resnode;

struct kdnode {
    size_t id;
    size_t id_orig;
    int dir;
    size_t left;
    size_t right;
};

struct kdtree {
    int ndim;
    size_t nnodes;
    size_t nallocated;
    kdnode* nodes;
    double* coords;
    double* min;
    double* max;
};

struct resnode {
    size_t id;
    double dist;
    resnode* next;
};

struct kdset {
    resnode* root;
    resnode* toread;
    size_t size;
    int beingread;
};

/**
 */
kdtree* kd_create(int ndim)
{
    kdtree* tree;
    int i;

    tree = malloc(sizeof(kdtree));
    tree->ndim = ndim;
    tree->nnodes = 0;
    tree->nallocated = 0;
    tree->nodes = NULL;
    tree->coords = NULL;
    tree->min = malloc(ndim * 2 * sizeof(double));
    tree->max = &tree->min[ndim];
    for (i = 0; i < ndim; ++i) {
        tree->min[i] = DBL_MAX;
        tree->max[i] = -DBL_MAX;
    }

    return tree;
}

/**
 */
void kd_destroy(kdtree* tree)
{
    if (tree->nallocated > 0) {
        free(tree->nodes);
        free(tree->coords);
    }
    free(tree->min);
    free(tree);
}

/**
 */
static int _kd_insertnode(kdtree* tree, size_t* id, const double* coords, size_t id_orig, int dir)
{
    int ndim = tree->ndim;
    int new_dir;
    kdnode* node;

    if (*id == SIZE_MAX) {
        node = &tree->nodes[tree->nnodes];
        memcpy(&tree->coords[tree->nnodes * ndim], coords, ndim * sizeof(double));
        node->id = tree->nnodes;
        node->id_orig = id_orig;
        node->dir = dir;
        node->left = SIZE_MAX;
        node->right = SIZE_MAX;
        *id = node->id;
        return 0;
    }
    node = &tree->nodes[*id];
    new_dir = (node->dir + 1) % ndim;

    return _kd_insertnode(tree, (coords[node->dir] < tree->coords[node->id * ndim + node->dir]) ? &node->left : &node->right, coords, id_orig, new_dir);
}

/**
 */
void kd_insertnode(kdtree* tree, const double* coords, size_t id_orig)
{
    int i;

    if (!isfinite(coords[0]))
        return;

    if (tree->nallocated == tree->nnodes) {
        if (tree->nallocated != 0)
            tree->nallocated *= 2;
        else
            tree->nallocated = NALLOCSTART;
        tree->nodes = realloc(tree->nodes, tree->nallocated * sizeof(kdnode));
        tree->coords = realloc(tree->coords, tree->nallocated * tree->ndim * sizeof(double));
        if (tree->nallocated == NALLOCSTART)
            tree->nodes[0].id = SIZE_MAX;
    }

    (void) _kd_insertnode(tree, &tree->nodes[0].id, coords, id_orig, 0);
    for (i = 0; i < tree->ndim; ++i) {
        if (coords[i] < tree->min[i])
            tree->min[i] = coords[i];
        if (coords[i] > tree->max[i])
            tree->max[i] = coords[i];
    }
    tree->nnodes++;
}

/**
 */
static void shuffle(size_t n, size_t ids[])
{
    size_t i;

    srand48(SEED);

    for (i = 0; i < n; ++i) {
        size_t ii = (size_t) ((double) n * drand48());
        size_t tmp = ids[i];

        ids[i] = ids[ii];
        ids[ii] = tmp;
    }
}

/**
 */
void kd_insertnodes(kdtree* tree, size_t n, double** src, int randomise)
{
    int nnodes0 = tree->nnodes;
    size_t* ids = NULL;
    double* coords;
    size_t i, j;

    if (n <= 0)
        return;

    if (tree->nallocated - tree->nnodes < n) {
        tree->nallocated += n;
        tree->nodes = realloc(tree->nodes, tree->nallocated * sizeof(kdnode));
        tree->coords = realloc(tree->coords, tree->nallocated * tree->ndim * sizeof(double));
        if (tree->nallocated == n)
            tree->nodes[0].id = SIZE_MAX;
    }

    if (randomise) {
        ids = malloc(n * sizeof(size_t));
        for (i = 0; i < n; ++i)
            ids[i] = i;
        shuffle(n, ids);
    }

    coords = malloc(tree->ndim * sizeof(double));

    for (i = 0; i < n; ++i) {
        int id = (randomise) ? ids[i] : i;

        for (j = 0; j < tree->ndim; ++j)
            coords[j] = src[j][id];
        if (!isfinite(coords[0]))
            continue;
        (void) _kd_insertnode(tree, &tree->nodes[0].id, coords, id + nnodes0, 0);
        tree->nnodes++;
        for (j = 0; j < tree->ndim; ++j) {
            if (coords[j] < tree->min[j])
                tree->min[j] = coords[j];
            if (coords[j] > tree->max[j])
                tree->max[j] = coords[j];
        }
    }

    if (randomise)
        free(ids);
    free(coords);
}

/**
 */
size_t kd_getsize(kdtree* tree)
{
    return tree->nnodes;
}

/**
 */
static double disttohyperrect(int ndim, const double* min, const double* max, const double* coords)
{
    double dist = 0.0;
    int i;

    for (i = 0; i < ndim; ++i) {
        if (coords[i] < min[i])
            dist += (min[i] - coords[i]) * (min[i] - coords[i]);
        else if (coords[i] > max[i])
            dist += (max[i] - coords[i]) * (max[i] - coords[i]);
    }

    return dist;
}

/**
 */
static void _kdset_insert(kdset* set, size_t id, double dist, int ordered)
{
    resnode* res = malloc(sizeof(resnode));

    res->id = id;
    res->dist = dist;

    if (set->root == NULL) {
        set->root = res;
        set->toread = NULL;
        set->size = 0;
        set->beingread = 0;
        res->next = NULL;
    } else if (ordered && set->root->dist > res->dist) {
        resnode* tmp = set->root;

        set->root = res;
        set->root->next = tmp;
    } else {
        resnode* now = set->root;

        if (ordered)
            while (now->next != NULL && now->next->dist < dist)
                now = now->next;
        res->next = now->next;
        now->next = res;
    }
}

/**
 */
static int _kd_findnodeswithinrange(const kdtree* tree, int id, const double* coords, double range, kdset* set, int ordered)
{
    int ndim = tree->ndim;
    kdnode* node;
    double* nodecoords;
    double dist, dx;
    int i, ret, added_res;

    if (id < 0)
        return 0;

    node = &tree->nodes[id];
    nodecoords = &tree->coords[node->id * ndim];

    for (i = 0, dist = 0.0; i < ndim; i++)
        dist += (nodecoords[i] - coords[i]) * (nodecoords[i] - coords[i]);

    added_res = 0;
    if (dist <= range * range) {
        _kdset_insert(set, id, dist, ordered);
        added_res = 1;
    }

    dx = coords[node->dir] - nodecoords[node->dir];
    ret = _kd_findnodeswithinrange(tree, dx <= 0.0 ? node->left : node->right, coords, range, set, ordered);
    if (fabs(dx) < range) {
        added_res += ret;
        ret = _kd_findnodeswithinrange(tree, dx <= 0.0 ? node->right : node->left, coords, range, set, ordered);
    }
    added_res += ret;

    return added_res;
}

/**
 */
kdset* kd_findnodeswithinrange(const kdtree* tree, const double* coords, double range, int ordered)
{
    int ret;
    kdset* rset;

    rset = malloc(sizeof(kdset));
    rset->root = NULL;
    rset->size = 0;

    ret = _kd_findnodeswithinrange(tree, 0, coords, range, rset, ordered);
    rset->size = ret;

    return rset;
}

/**
 */
static void _kd_findnearestnode(const kdtree* tree, const size_t nodeid, const double* coords, size_t* result, double* resdist, double* minmax)
{
    int ndim = tree->ndim;
    kdnode* node = &tree->nodes[nodeid];
    double* nodecoords = &tree->coords[nodeid * ndim];
    int dir = node->dir;
    int left = (coords[dir] - nodecoords[dir] <= 0.0);
    size_t nearer_subtree, farther_subtree;
    double* nearer_hyperrect_coord;
    double* farther_hyperrect_coord;
    double dist;
    int i;

    if (left) {
        nearer_subtree = node->left;
        farther_subtree = node->right;
        nearer_hyperrect_coord = &minmax[ndim + dir];
        farther_hyperrect_coord = &minmax[dir];
    } else {
        nearer_subtree = node->right;
        farther_subtree = node->left;
        nearer_hyperrect_coord = &minmax[dir];
        farther_hyperrect_coord = &minmax[ndim + dir];
    }

    if (nearer_subtree != SIZE_MAX) {
        /*
         * get the hyperrect of the nearer subtree 
         */
        double tmp = nearer_hyperrect_coord[0];

        nearer_hyperrect_coord[0] = nodecoords[dir];
        /*
         * recurse down into nearer subtree 
         */
        _kd_findnearestnode(tree, nearer_subtree, coords, result, resdist, minmax);
        /*
         * undo the slice 
         */
        nearer_hyperrect_coord[0] = tmp;
    }

    /*
     * check the distance of the point at the current node, compare it with
     * the best so far 
     */
    for (i = 0, dist = 0.0; i < ndim; ++i)
        dist += (nodecoords[i] - coords[i]) * (nodecoords[i] - coords[i]);
    if (dist <= *resdist) {
        *result = nodeid;
        *resdist = dist;
    }

    if (farther_subtree != SIZE_MAX) {
        /*
         * get the hyperrect of the farther subtree 
         */
        double tmp = farther_hyperrect_coord[0];

        farther_hyperrect_coord[0] = nodecoords[dir];
        /*
         * check if we have to recurse down by calculating the closest point
         * of the hyperrect and see if it's closer than our minimum distance
         * in result_dist_sq. 
         */
        if (disttohyperrect(ndim, minmax, &minmax[ndim], coords) < *resdist)
            /*
             * recurse down into farther subtree 
             */
            _kd_findnearestnode(tree, farther_subtree, coords, result, resdist, minmax);
        /*
         * undo the slice 
         */
        farther_hyperrect_coord[0] = tmp;
    }
}

/**
 */
size_t kd_findnearestnode(const kdtree* tree, const double* coords)
{
    int ndim = tree->ndim;
    double* minmax = malloc(ndim * 2 * sizeof(double));
    size_t result;
    double dist;
    int i;

    for (i = 0; i < ndim * 2; ++i)
        minmax[i] = tree->min[i];

    for (i = 0, dist = 0.0; i < ndim; ++i)
        dist += (tree->coords[i] - coords[i]) * (tree->coords[i] - coords[i]);

    /*
     * search for the nearest neighbour recursively 
     */
    _kd_findnearestnode(tree, 0, coords, &result, &dist, minmax);

    free(minmax);

    return result;
}

/**
 */
static void _clear_results(kdset* rset)
{
    resnode* root = rset->root;

    while (root != NULL) {
        resnode* now = root;

        root = root->next;
        free(now);
    }
    rset->root = NULL;
    rset->size = 0;
}

/**
 */
void kdset_free(kdset* set)
{
    _clear_results(set);
    free(set);
}

/**
 */
size_t kdset_getsize(const kdset* set)
{
    return set->size;
}

/**
 */
size_t kdset_read(kdset* set, double* dist)
{
    size_t id;

    if (set == NULL || set->root == NULL) {
        *dist = NAN;
        return SIZE_MAX;
    }

    if (set->beingread == 0) {
        set->beingread = 1;
        set->toread = set->root;
    }

    if (set->toread == NULL) {
        *dist = NAN;
        return SIZE_MAX;
    }

    *dist = sqrt(set->toread->dist);
    id = set->toread->id;
    set->toread = set->toread->next;

    return id;
}

/**
 */
double* kd_getnodecoords(const kdtree* tree, size_t id)
{
    return &tree->coords[id * tree->ndim];
}

/**
 */
size_t kd_getnodeorigid(const kdtree* tree, size_t id)
{

    return (int) tree->nodes[id].id_orig;
}

/* get boundary rectangle
 */
double* kd_getminmax(const kdtree* tree)
{
    return tree->min;
}

#if defined(STANDALONE)

#include <stdarg.h>
#include <errno.h>
#include "nan.h"

#define BUFSIZE 1024
#define NDIMMAX 256
#define NALLOCSTART 1024
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
static kdtree* kd_readfile(char* fname, int ndimin, int dorandomise)
{
    int ndim = (isll) ? 3 : ndim;
    FILE* f = NULL;
    kdtree* tree = NULL;
    char buf[BUFSIZE];
    char seps[] = " ,;\t";
    char* token;
    double** coords;
    double point[NDIMMAX];
    size_t nallocated;
    size_t npoints;
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

    coords = malloc(ndim * sizeof(double*));
    nallocated = NALLOCSTART;
    for (i = 0; i < ndim; ++i)
        coords[i] = malloc(nallocated * sizeof(double));

    npoints = 0;
    while (fgets(buf, BUFSIZE, f) != NULL) {
        if (buf[0] == '#')
            continue;

        for (i = 0; i < ndimin; ++i) {
            if ((token = strtok((i == 0) ? buf : NULL, seps)) == NULL)
                break;
            if (!str2double(token, &point[i]))
                break;
        }
        if (i < ndimin)
            continue;

        if (isll) {
            double lon = point[0] * INVRAD;
            double lat = point[1] * INVRAD;
            double coslat = cos(lat);

            point[0] = REARTH * sin(lon) * coslat;
            point[1] = REARTH * cos(lon) * coslat;
            point[2] = REARTH * sin(lat);
        }

        if (npoints == nallocated) {
            nallocated *= 2;
            for (i = 0; i < ndim; ++i)
                coords[i] = realloc(coords[i], nallocated * sizeof(double));
        }

        for (i = 0; i < ndim; ++i)
            coords[i][npoints] = point[i];
        npoints++;
    }
    if (f != stdin)
        if (fclose(f) != 0)
            quit("%s: %s", fname, strerror(errno));

    tree = kd_create(ndim);
    kd_insertnodes(tree, npoints, coords, dorandomise);

    for (i = 0; i < ndim; ++i)
        free(coords[i]);
    free(coords);

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
    printf("  Usage: kd -i <file> -n <ndim> -p <coords> -r <range> [-s] [-g] [-o]\n");
    printf("         kd -h\n");
    printf("  Options:\n");
    printf("    -i <file>   -- a text file, with each row containing <ndim> point\n");
    printf("                   coordinates\n");
    printf("    -g          -- assume coordinates to be lon/lat in degrees and the distance is in\n");
    printf("                   kilometres (ndim = 2)\n");
    printf("    -h          -- describe this program\n");
    printf("    -n <ndim>   -- number of dimensions\n");
    printf("    -o          -- do not shuffle input points\n");
    printf("    -p <coords> -- <ndim> coordinates of the point to search around\n");
    printf("    -r <range>  -- radius to search within\n");
    printf("    -s          -- sort results\n");
    printf("                   <range> in km\n");

    exit(0);
}

/**
 */
static void parse_commandline(int argc, char* argv[], char** fname, int* ndim, double* coords, double* range, int* dosort, int* dorandomise)
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
        case 'g':
            i++;
            isll = 1;
            if (*ndim <= 0)
                *ndim = 2;
            break;
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
            for (ii = 0; i < argc && argv[i][0] != '-'; ++ii, ++i)
                if (!str2double(argv[i], &coords[ii]))
                    usage();
            break;
        case 'o':
            i++;
            *dorandomise = 0;
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
        case 'h':
            description();
            exit(0);
        default:
            usage();
        }
    }

    if (*ndim <= 0)
        quit("\"ndim\" must be defined");
}

/**
 */
int main(int argc, char* argv[])
{
    char* fname = NULL;
    int ndim = -1;
    double coords[NDIMMAX];
    double range = NaN;
    int dosort = 0;
    int dorandomise = 1;
    kdtree* tree = NULL;
    kdset* set = NULL;
    double dist;
    size_t id, id_orig;
    int i;

    parse_commandline(argc, argv, &fname, &ndim, coords, &range, &dosort, &dorandomise);

    tree = kd_readfile(fname, ndim, dorandomise);
    if (kd_getsize(tree) == 0) {
        printf("  nothing to do, exiting\n");
        return 0;
    }

    if (isll) {
        double lon = coords[0] * INVRAD;
        double lat = coords[1] * INVRAD;

        coords[0] = REARTH * sin(lon) * cos(lat);
        coords[1] = REARTH * cos(lon) * cos(lat);
        coords[2] = REARTH * sin(lat);
    }

    /*
     * first find the nearest node
     */
    printf("  searching for the nearest node:");
    fflush(stdout);
    id = kd_findnearestnode(tree, coords);
    id_orig = kd_getnodeorigid(tree, id);
    printf("\n");
    printf("       X        Y       ID      DIST\n");
    {
        double* rescoords = kd_getnodecoords(tree, id);
        double rescoordsout[NDIMMAX];

        if (isll) {
            double lat = asin(rescoords[2] / REARTH);
            double lon = asin(rescoords[0] / REARTH / cos(lat));

            dist = sqrt((rescoords[0] - coords[0]) * (rescoords[0] - coords[0]) + (rescoords[1] - coords[1]) * (rescoords[1] - coords[1]) + (rescoords[2] - coords[2]) * (rescoords[2] - coords[2]));
            rescoordsout[0] = lon / INVRAD;
            rescoordsout[1] = lat / INVRAD;
        } else {
            for (i = 0; i < ndim; ++i)
                dist += (rescoords[i] - coords[i]) * (rescoords[i] - coords[i]);
            for (i = 0; i < ndim; ++i)
                rescoordsout[i] = rescoords[i];
            dist = sqrt(dist);
        }
        printf("  ");
        for (i = 0; i < ndim; ++i)
            printf("%8.3f ", rescoordsout[i]);
        printf("%8zu ", id_orig);
        printf("%8g\n", dist);
    }
    fflush(stdout);

    /*
     * first nodes within range
     */
    printf("  searching for the nodes within range:");
    fflush(stdout);
    set = kd_findnodeswithinrange(tree, coords, range, dosort);
    printf("\n");
    printf("    %zu points found:\n", kdset_getsize(set));
    if (kdset_getsize(set) > 0)
        printf("       X        Y       ID      DIST\n");
    for (; (id = kdset_read(set, &dist)) < SIZE_MAX;) {
        double* rescoords = kd_getnodecoords(tree, id);

        id_orig = kd_getnodeorigid(tree, id);

        if (isll) {
            double lat = asin(rescoords[2] / REARTH);
            double lon = asin(rescoords[0] / REARTH / cos(lat));

            rescoords[0] = lon / INVRAD;
            rescoords[1] = lat / INVRAD;
        }
        printf("  ");
        for (i = 0; i < ndim; ++i)
            printf("%8.3f ", rescoords[i]);
        printf("%8zu ", id_orig);
        printf("%8g\n", dist);
    }

    kdset_free(set);
    kd_destroy(tree);

    return 0;
}
#endif                          /* STANDALONE */

/*
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
