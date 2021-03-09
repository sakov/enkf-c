/******************************************************************************
 *
 * File:        kdtree.c        
 *
 * Created:     23/03/2016
 *
 * Author:      Pavel Sakov
 *
 * Description: KD-tree code.  Initially derived from the code by John
 *              Tsiombikas (see the tail of the file). This implementation is
 *              driven by the needs of large scale geophysical data
 *              assimilation systems. The main changes concern (1) using
 *              continuous block of memory and (2) making it possible to
 *              pre-allocate memory externally to permit using shared memory.
 *
 * Revisions:   23/11/2018 PS:
 *              - Replaced individual allocations of "resnode" in the linked
 *                list structure in _kdset_insert() by allocations of blocks of
 *                KDSET_BLOCKSIZE resnodes.
 *              - Streamlined counting the size of the results set by simply
 *                adding "set->size++" to _kdset_insert().
 *              6/12/2019 PS:
 *              - The field "id_orig" became "data". In kd_insertnode() it is
 *                treated as generic data; in kd_insertnodes() it can be either
 *                generic data or (if NULL) the sequential number of the node.
 *              27/3/2020 PS:
 *              - Added  kd_getstoragesize(), kd_setstorage(), and
 *                kd_syncsize().
 *              6/3/2020 PS:
 *              - A number of changes for range search targeting large systems;
 *                in particular, to minimise the number of memory allocations
 *                for multiple repeated searches. The results are now stored
 *                within the tree to re-use the allocated memory and are only
 *                deallocated in the destructor.
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <stdarg.h>
#include <errno.h>
#include <assert.h>
#include "kdtree.h"

#define NALLOCSTART 1024
#define NDIMMAX 3
#define RES_INC 4096

typedef struct {
    size_t id;
    size_t data;                /* either set explicitly when adding node(s)
                                 * or set to the sequential number of the
                                 * node added */
    int dir;
    size_t left;
    size_t right;
} kdnode;

struct kdtree {
    char* name;
    size_t ndim;
    size_t nnodes;
    size_t nallocated;
    kdnode* nodes;
    double* coords;
    double* min;
    double* max;
    int external_storage;

    size_t nresults;
    size_t nresults_allocated;
    kdresult* results;
};

/**
 */
static void quit(char* format, ...)
{
    va_list args;

    fflush(stdout);

    fprintf(stderr, "\n\n  ERROR: kdtree: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n\n");
    fflush(NULL);
    exit(1);
}

/**
 */
kdtree* kd_create(char* name, size_t ndim)
{
    kdtree* tree;
    int i;

    if (ndim > NDIMMAX)
        quit("ndim = %d, NDIMMAX = %d", ndim, NDIMMAX);
    tree = malloc(sizeof(kdtree));
    tree->name = strdup(name);
    tree->ndim = ndim;
    tree->nnodes = 0;
    tree->nallocated = 0;
    tree->min = malloc(ndim * 2 * sizeof(double));
    tree->nodes = (void*) tree->min;
    tree->coords = (void*) tree->min;
    tree->max = &tree->min[ndim];
    tree->external_storage = 0;

    for (i = 0; i < ndim; ++i) {
        tree->min[i] = DBL_MAX;
        tree->max[i] = -DBL_MAX;
    }

    tree->nresults = 0;
    tree->nresults_allocated = 0;
    tree->results = NULL;

    return tree;
}

/**
 */
void kd_destroy(kdtree* tree)
{
    if (tree->nallocated > 0 && !tree->external_storage)
        free(tree->nodes);
    free(tree->name);
    if (tree->nresults_allocated > 0)
        free(tree->results);
    free(tree);
}

/**
 */
static int _kd_insertnode(kdtree* tree, size_t* id, const double* coords, size_t data, int dir)
{
    int ndim = tree->ndim;
    int new_dir;
    kdnode* node;

    if (*id == SIZE_MAX) {
        node = &tree->nodes[tree->nnodes];
        memcpy(&tree->coords[tree->nnodes * ndim], coords, ndim * sizeof(double));
        node->id = tree->nnodes;
        node->data = data;
        node->dir = dir;
        node->left = SIZE_MAX;
        node->right = SIZE_MAX;
        *id = node->id;
        return 0;
    }
    node = &tree->nodes[*id];
    new_dir = (node->dir + 1) % ndim;

    return _kd_insertnode(tree, (coords[node->dir] < tree->coords[node->id * ndim + node->dir]) ? &node->left : &node->right, coords, data, new_dir);
}

static void _kd_changealloc(kdtree* tree, size_t nallocated)
{
    int ndim = tree->ndim;

    assert(nallocated >= tree->nnodes);
    if (nallocated > tree->nallocated) {
        double* coords_prev;
        double* min_prev;

        tree->nodes = realloc(tree->nodes, nallocated * sizeof(kdnode) + nallocated * ndim * sizeof(double) + ndim * 2 * sizeof(double));
        tree->coords = (double*) &tree->nodes[nallocated];
        coords_prev = (double*) &tree->nodes[tree->nallocated];
        tree->min = &tree->coords[nallocated * ndim];
        min_prev = &coords_prev[tree->nallocated * ndim];
        memmove(tree->min, min_prev, ndim * 2 * sizeof(double));
        memmove(tree->coords, coords_prev, tree->nnodes * ndim * sizeof(double));
        tree->max = &tree->min[ndim];
        if (tree->nnodes == 0)
            tree->nodes[0].id = SIZE_MAX;
    } else if (nallocated < tree->nallocated) {
        tree->coords = (double*) &tree->nodes[nallocated];
        memmove(tree->coords, &tree->nodes[tree->nallocated], tree->nnodes * ndim * sizeof(double));
        tree->min = &tree->coords[nallocated * ndim];
        memmove(tree->min, &tree->coords[tree->nallocated * ndim], ndim * 2 * sizeof(double));
        tree->nodes = realloc(tree->nodes, nallocated * sizeof(kdnode) + nallocated * ndim * sizeof(double) + ndim * 2 * sizeof(double));
        tree->coords = (double*) &tree->nodes[nallocated];
        tree->min = &tree->coords[nallocated * ndim];
        tree->max = &tree->min[ndim];
    }
    tree->nallocated = nallocated;
}

/**
 */
void kd_insertnode(kdtree* tree, const double* coords, size_t data)
{
    size_t ndim = tree->ndim;
    size_t i;

    assert(tree->nallocated >= tree->nnodes);

    if (!isfinite(coords[0]))
        return;

    if (tree->nallocated == tree->nnodes)
        _kd_changealloc(tree, (tree->nallocated != 0) ? tree->nallocated * 2 : NALLOCSTART);

    (void) _kd_insertnode(tree, &tree->nodes[0].id, coords, data, 0);
    for (i = 0; i < ndim; ++i) {
        if (coords[i] < tree->min[i])
            tree->min[i] = coords[i];
        if (coords[i] > tree->max[i])
            tree->max[i] = coords[i];
    }
    tree->nnodes++;
}

/**
 */
static void randomise_rand48(void)
{
    char fname[] = "/dev/urandom";
    FILE* f = NULL;
    size_t status;

    if (seed_rand48 != 0)
        return;

    f = fopen(fname, "r");
    if (f == NULL) {
        int errno_saved = errno;

        quit("randomise_rand48(): could not open \"%s\": %s", fname, strerror(errno_saved));
    }

    status = fread(&seed_rand48, sizeof(seed_rand48), 1, f);
    if (status != 1) {
        int errno_saved = errno;

        quit("randomise_rand48(): could not read from \"%s\": %s", fname, strerror(errno_saved));
    }
    fclose(f);

    srand(seed_rand48);
}

/**
 */
static void shuffle(size_t n, size_t ids[])
{
    size_t i;

    /*
     * (PS) this initialisation is not really necessary, I think, but let it
     * stay for now
     */
    randomise_rand48();

    for (i = 0; i < n; ++i) {
        size_t ii = (size_t) ((double) n * drand48());
        size_t tmp = ids[i];

        ids[i] = ids[ii];
        ids[ii] = tmp;
    }
}

/** Insert nodes into the tree.
 * @param tree - Kd-tree
 * @param n - number of inserted nodes
 * @param src - [ndim][n] arrays of coordinates of inserted nodes
 * @param data - [n] array of data or NULL
 * @param mask - [n] array of mask values (node [i] is ignored if mask[i] = 0)
 * @param randomise - whether to insert nodes in random order. Should normally
 *                    be set to 1 for getting a balanced tree.
 */
void kd_insertnodes(kdtree* tree, size_t n, double** src, size_t* data, int* mask, int randomise)
{
    size_t nnodes_prev = tree->nnodes;
    size_t* ids = NULL;
    double* coords;
    size_t i, j, ngood;

    assert(tree->nallocated >= tree->nnodes);

    if (n <= 0)
        return;

    if (tree->nallocated - tree->nnodes < n)
        _kd_changealloc(tree, tree->nnodes + n);

    if (randomise || mask != NULL) {
        ids = malloc(n * sizeof(size_t));
        for (i = 0, ngood = 0; i < n; ++i) {
            if (mask != NULL && mask[i] == 0)
                continue;
            ids[ngood] = i;
            ngood++;
        }
        if (randomise)
            shuffle(ngood, ids);
    } else
        ngood = n;

    coords = malloc(tree->ndim * sizeof(double));

    for (i = 0; i < ngood; ++i) {
        int id = (ids != NULL) ? ids[i] : i;

        for (j = 0; j < tree->ndim; ++j)
            coords[j] = src[j][id];
        if (!isfinite(coords[0]))
            continue;
        (void) _kd_insertnode(tree, &tree->nodes[0].id, coords, (data != NULL) ? data[id] : nnodes_prev + id, 0);
        tree->nnodes++;
        for (j = 0; j < tree->ndim; ++j) {
            if (coords[j] < tree->min[j])
                tree->min[j] = coords[j];
            if (coords[j] > tree->max[j])
                tree->max[j] = coords[j];
        }
    }

    if (ids != NULL)
        free(ids);
    free(coords);
}

/** Allocate space for n nodes (for an empty tree only).
 * @param tree     Kd-tree
 * @param n        Number of nodes to allocate for
 * @param storage  The external storage to use (must have the size to hold n
 *                 nodes).
 * @param ismaster Flag whether to initialise min/max
 */
void kd_setstorage(kdtree* tree, size_t n, void* storage, int ismaster)
{
    size_t ndim = tree->ndim;

    assert(tree->nnodes == 0);
    assert(n > 0);

    if (storage == NULL) {
        _kd_changealloc(tree, n);
        return;
    }

    free(tree->nodes);

    tree->nodes = storage;
    tree->coords = (double*) &tree->nodes[n];
    tree->min = &tree->coords[ndim * n];
    tree->max = &tree->min[ndim];
    if (ismaster) {
        size_t i;

        tree->nodes[0].id = SIZE_MAX;
        for (i = 0; i < ndim; ++i) {
            tree->min[i] = DBL_MAX;
            tree->max[i] = -DBL_MAX;
        }
    }

    tree->nallocated = n;
    tree->external_storage = 1;
}

/**
 */
void kd_syncsize(kdtree* tree)
{
    assert(tree->nnodes == 0);
    tree->nnodes = tree->nallocated;
}

/**
 */
void kd_finalise(kdtree* tree)
{
    if (tree->nallocated == tree->nnodes)
        return;

    memmove(&tree->nodes[tree->nnodes], tree->coords, tree->nnodes * tree->ndim * sizeof(double));
    tree->nodes = realloc(tree->nodes, tree->nnodes * sizeof(kdnode) + tree->nnodes * tree->ndim * sizeof(double));
    tree->coords = (double*) &tree->nodes[tree->nnodes];
    tree->nallocated = tree->nnodes;
}

/**
 */
size_t kd_getsize(kdtree* tree)
{
    return tree->nnodes;
}

/**
 */
size_t kd_getnalloc(kdtree* tree)
{
    return tree->nallocated;
}

/**
 */
size_t kd_getnodesize(kdtree* tree)
{
    return sizeof(kdnode) + tree->ndim * sizeof(double);
}

/**
 */
char* kd_getname(const kdtree* tree)
{
    return tree->name;
}

/** get the memory size necessary for the tree for the current tree
 *  (nnodes = 0) or for the tree with `nnodes' nodes
 */
size_t kd_getstoragesize(const kdtree* tree, size_t nnodes)
{
    if (nnodes == 0)
        nnodes = tree->nnodes;

    return nnodes * sizeof(kdnode) + nnodes * tree->ndim * sizeof(double) + tree->ndim * 2 * sizeof(double);
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
static void _insertresult(kdtree* tree, size_t id, double distsq)
{
    kdresult* res;

    if (tree->nresults >= tree->nresults_allocated) {
        tree->nresults_allocated += RES_INC;
        tree->results = realloc(tree->results, tree->nresults_allocated * sizeof(kdresult));
    }

    res = &tree->results[tree->nresults];
    res->id = id;
    res->distsq = distsq;
    tree->nresults++;
}

/**
 */
static void _kd_findnodeswithinrange(kdtree* tree, int id, const double* coords, double range)
{
    int ndim = tree->ndim;
    kdnode* node;
    double* nodecoords;
    double distsq, dx;
    int i;

    if (id < 0)
        return;

    node = &tree->nodes[id];
    nodecoords = &tree->coords[node->id * ndim];

    for (i = 0, distsq = 0.0; i < ndim; i++)
        distsq += (nodecoords[i] - coords[i]) * (nodecoords[i] - coords[i]);

    if (distsq <= range * range)
        _insertresult(tree, id, distsq);

    dx = coords[node->dir] - nodecoords[node->dir];
    _kd_findnodeswithinrange(tree, dx <= 0.0 ? node->left : node->right, coords, range);
    if (fabs(dx) < range)
        _kd_findnodeswithinrange(tree, dx <= 0.0 ? node->right : node->left, coords, range);
}

/**
 */
static int cmp_res(const void* p1, const void* p2)
{
    kdresult* n1 = (kdresult *) p1;
    kdresult* n2 = (kdresult *) p2;

    if (n1->distsq > n2->distsq)
        return 1;
    else if (n1->distsq < n2->distsq)
        return -1;
    return 0;
}

/**
 */
void kd_findnodeswithinrange(kdtree* tree, const double* coords, double range, int ordered, size_t* nresults, kdresult ** results)
{
    tree->nresults = 0;
    _kd_findnodeswithinrange(tree, 0, coords, range);

    if (ordered)
        qsort(tree->results, tree->nresults, sizeof(kdresult), cmp_res);

    *nresults = tree->nresults;
    *results = tree->results;
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

/** Find the nearest node to the specified point
 * @param tree - Kd-tree
 * @param coords - [ndim] array of point coordinates
 * @return - ID of the nearest node
 */
size_t kd_findnearestnode(const kdtree* tree, const double* coords)
{
    int ndim = tree->ndim;
    double minmax[NDIMMAX * 2];
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

    return result;
}

/** Get coordinates of the node.
 * @param tree - Kd-tree
 * @param id - ID of the node
 * @return - [ndim] array of node coordinates
 */
double* kd_getnodecoords(const kdtree* tree, size_t id)
{
    if (id >= tree->nnodes)
        quit("id = %zu >= tree size = %zu", id, tree->nnodes);

    return &tree->coords[id * tree->ndim];
}

/**
 */
size_t kd_getnodedata(const kdtree* tree, size_t id)
{
    if (id >= tree->nnodes)
        quit("id = %zu >= tree size = %zu", id, tree->nnodes);

    return (int) tree->nodes[id].data;
}

/** get boundary rectangle
 */
double* kd_getminmax(const kdtree* tree)
{
    return tree->min;
}

#if defined(STANDALONE)

#define BUFSIZE 1024
#define NALLOCSTART 1024
#define NSRC_INC 10

#define INVRAD (3.14159265358979323846 / 180.0)
#define REARTH 6371.0

static int isll = 0;
long int seed_rand48 = 0;

/**
 */
static int str2double(char* token, double* value)
{
    char* end = NULL;

    if (token == NULL) {
        *value = NAN;
        return 0;
    }

    *value = strtod(token, &end);

    if (end == token) {
        *value = NAN;
        return 0;
    }

    return 1;
}

/**
 */
static void kd_readfile(kdtree* tree, char* fname, int ndimin)
{
    int ndim = (isll) ? 3 : ndimin;
    FILE* f = NULL;
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

            point[0] = REARTH * cos(lon) * coslat;
            point[1] = REARTH * sin(lon) * coslat;
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

    kd_insertnodes(tree, npoints, coords, NULL, NULL, 1);

    for (i = 0; i < ndim; ++i)
        free(coords[i]);
    free(coords);
}

/**
 */
static void description(void)
{
    printf("  kd: reads coordinates from text files; searches either for points\n");
    printf("      within specified radius or for the nearest point to the specified\n");
    printf("      location; prints results (point coordinates and ids) to the standard\n");
    printf("      output. The data can be either in cartesian (dimension <ndim>) or\n");
    printf("      geographic (lon,lat) coordinates.\n");
}

/**
 */
static void usage(void)
{
    printf("  Usage: kd -i <file> [<file> ...] {-n <ndim> | -g} -p <coords> [-r <range>] [-s]\n");
    printf("         kd -h\n");
    printf("  Options:\n");
    printf("    -i <file>   -- a text file, with each row containing <ndim> points\n");
    printf("                   coordinates\n");
    printf("    -g          -- assume coordinates to be lon/lat in degrees and distance to be in\n");
    printf("                   kilometres\n");
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
static void parse_commandline(int argc, char* argv[], int* nsrc, char*** srcs, int* ndim, double* coords, double* range, int* dosort)
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
            while (i < argc && argv[i][0] != '-') {
                if (*nsrc % NSRC_INC == 0)
                    *srcs = realloc(*srcs, (*nsrc + NSRC_INC) * sizeof(void*));
                (*srcs)[*nsrc] = strdup(argv[i]);
                (*nsrc)++;
                i++;
            }
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
            if (*ndim <= 0)
                quit("-n needs to be entered before -p");
            for (ii = 0; ii < *ndim; ++ii, ++i)
                if (!str2double(argv[i], &coords[ii]))
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
static void print_location(int ndim, double coords[])
{
    int i;

    printf("  location = (");
    for (i = 0; i < ndim; ++i)
        printf("%s%.6g", (i > 0) ? ", " : "", coords[i]);
    printf(")\n");
}

/**
 */
int main(int argc, char* argv[])
{
    int nsrc = 0;
    char** srcs = NULL;
    int ndim = -1;
    double coords[NDIMMAX];
    double range = NAN;
    int dosort = 0;

    kdtree* tree = NULL;
    double dist;
    int i;

    parse_commandline(argc, argv, &nsrc, &srcs, &ndim, coords, &range, &dosort);

    printf("  srcs = ");
    for (i = 0; i < nsrc; ++i)
        printf("%s ", srcs[i]);
    printf("\n");

    tree = kd_create("kdtree", (isll) ? 3 : ndim);
    for (i = 0; i < nsrc; ++i)
        kd_readfile(tree, srcs[i], ndim);
    printf("    %zu points\n", kd_getsize(tree));
    if (kd_getsize(tree) == 0) {
        printf("  nothing to do, exiting\n");
        return 0;
    }

    print_location(ndim, coords);
    if (isll) {
        double lon = coords[0] * INVRAD;
        double lat = coords[1] * INVRAD;

        coords[0] = REARTH * cos(lon) * cos(lat);
        coords[1] = REARTH * sin(lon) * cos(lat);
        coords[2] = REARTH * sin(lat);
    }

    if (isnan(range)) {
        size_t id;

        /*
         * find the nearest node
         */
        printf("  the nearest node:\n");
        fflush(stdout);
        id = kd_findnearestnode(tree, coords);
        printf("       X        Y       ID      DIST\n");
        {
            double* rescoords = kd_getnodecoords(tree, id);
            double rescoordsout[NDIMMAX];

            if (isll) {
                double lat = asin(rescoords[2] / REARTH);
                double lon = atan2(rescoords[1] / REARTH, rescoords[0] / REARTH);

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
            printf("%8zu ", kd_getnodedata(tree, id));
            printf("%8.3f\n", dist);
        }
        fflush(stdout);
    } else {
        size_t nresults;
        kdresult* results;
        size_t i;

        /*
         * find nodes within range
         */
        printf("  nodes within r < %.3g%s:\n", range, (isll) ? "km" : "");
        fflush(stdout);
        kd_findnodeswithinrange(tree, coords, range, dosort, &nresults, &results);
        printf("    %zu points found:\n", nresults);
        if (nresults > 0)
            printf("       X        Y       ID      DIST\n");
        for (i = 0; i < nresults; ++i) {
            size_t id = results[i].id;
            double dist = sqrt(results[i].dist);
            double* coords = kd_getnodecoords(tree, id);
            int idim;

            if (isll) {
                double lat = asin(coords[2] / REARTH);
                double lon = atan2(coords[1] / REARTH, coords[0] / REARTH);

                coords[0] = lon / INVRAD;
                coords[1] = lat / INVRAD;
            }

            printf("  ");
            for (idim = 0; idim < ndim; ++idim)
                printf("%8.3f ", coords[idim]);
            printf("%8zu ", kd_getnodedata(tree, id));
            printf("%8.3f\n", dist);
        }
    }

    for (i = 0; i < nsrc; ++i)
        free(srcs[i]);
    free(srcs);

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
