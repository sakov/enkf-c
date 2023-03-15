/******************************************************************************
 *
 * File:        kdtree.h        
 *
 * Created:     23/03/2016
 *
 * Author:      Pavel Sakov
 *
 * Description: KD-tree code.  Initially derived from the code by John
 *              Tsiombikas (see the tail of the file). The main changes concern
 *              (1) using continuous block of memory and (2) making it possible
 *              to pre-allocate memory externally to permit using shared memory.
 *
 * Revisions:   
 *
 *****************************************************************************/

#if !defined(_KDTREE_H_)

extern long int seed_rand48;

struct kdtree;
typedef struct kdtree kdtree;

/* It is possible to set the quit procedure. By default, the internal procedure
 * is used.
 */
typedef void (*kd_quit_fn) (char* format, ...);
void kd_set_quitfn(kd_quit_fn quit_fn);

typedef struct {
    size_t id;
    double distsq;
} kdresult;

/*
 * basic procedures
 */
kdtree* kd_create(char* name, size_t ndim);
void kd_destroy(kdtree* tree);
void kd_insertnode(kdtree* tree, const double* coords, size_t data);
void kd_insertnodes(kdtree* tree, size_t n, double** src, size_t* data, int* mask, int randomise);
void kd_finalise(kdtree* tree);

/*
 * memory pre-allocation
 */
size_t kd_getstoragesize(const kdtree* tree, size_t nnodes);
void kd_setstorage(kdtree* tree, size_t n, void* storage, int ismaster);
void kd_syncsize(kdtree* tree);

/*
 * tree information
 */
char* kd_getname(const kdtree* tree);
size_t kd_getsize(kdtree* tree);
size_t kd_getnalloc(kdtree* tree);
size_t kd_getndim(kdtree* tree);
double* kd_getminmax(const kdtree* tree);
double* kd_getnodecoords(const kdtree* tree, size_t id);
size_t kd_getnodedata(const kdtree* tree, size_t id);

/*
 * searches
 */
void kd_findnodeswithinrange(kdtree* tree, const double* coords, double range, int ordered, size_t* n, kdresult** results);
size_t kd_findnearestnode(const kdtree* tree, const double* coords);

#define _KDTREE_H_
#endif                          /* _KDTREE_H_ */

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
