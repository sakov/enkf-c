/******************************************************************************
 *
 * File:        kdtree.h        
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

#if !defined(_KDTREE_H_)

extern long int seed_rand48;

struct kdtree;
typedef struct kdtree kdtree;

struct kdnode;
typedef struct kdnode kdnode;

struct kdset;
typedef struct kdset kdset;

/* create a kd-tree for "k"-dimensional data
 */
kdtree* kd_create(int ndim);

/* free the kdtree
 */
void kd_destroy(kdtree* tree);

/* insert a node, specifying its position, and optional data
 */
void kd_insertnode(kdtree* tree, const double* coords, size_t id_extern);

/* insert an array of nodes
 */
void kd_insertnodes(kdtree* tree, size_t n, double** src, size_t* ids_extern, int* mask, int randomise);

/* get the number of tree nodes
 */
size_t kd_getsize(kdtree* tree);

/* find any nearest nodes from the specified point within a range
 */
kdset* kd_findnodeswithinrange(const kdtree* tree, const double* coords, double range, int ordered);

/* find the nearest node
 */
size_t kd_findnearestnode(const kdtree* tree, const double* coords);

/* get position of a node
 */
double* kd_getnodecoords(const kdtree* tree, size_t id);

/* get the external ID of the node (it is different to the current ID if the
 * input was shuffled)
 */
size_t kd_getnodedata(const kdtree* tree, size_t id);

/* get boundary rectangle
 */
double* kd_getminmax(const kdtree* tree);

/* read node id of the current result (SIZE_MAX if no more results are
 * available; advance the result set iterator)
 */
size_t kdset_readnext(kdset* set, double* dist);

/* get the size of the result set
 */
size_t kdset_getsize(const kdset* set);

/* free a result set
 */
void kdset_free(kdset* set);

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
