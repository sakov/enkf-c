/******************************************************************************
 *
 * File:        kdtree.h        
 *
 * Created:     12/2012
 *
 * Author:      John Tsiombikas <nuclear@siggraph.org>
 *
 * Description:
 *
 * Revisions:   April 2012 - modified by Pavel Sakov
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

#if !defined(_KDTREE_H_)

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
int kd_insert(kdtree* tree, const double* pos);

/* find any nearest nodes from the specified point within a range
 */
kdset* kd_nearest_range(kdtree* tree, const double* pos, double range, int ordered);

/* get position of a node
 */
double* kd_getpos(kdtree* tree, int id);

/* free a result set
 */
void kd_res_free(kdset* set);

/* advance the result set iterator
 */
int kd_res_next(kdset* set);

/* get the current node ID in the result set
 */
int kd_res_getid(kdset* set);

#define _KDTREE_H_
#endif                          /* _KDTREE_H_ */
