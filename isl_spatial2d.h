#ifndef ISL_INCLUDE_SPATIAL2D_H_
#define ISL_INCLUDE_SPATIAL2D_H_

/* isl_spatial2d - v0.1 - public domain spatial hashing in 2d space
 * no warranty implied; use at your own risk
 *
 * Do this:
 *   #define ISL_SPATIAL2D_IMPLEMENTATION
 * before you include this file in *one* C or C++ file to create the implementation.
 * Internally uses stb_ds.h for dynamic arrays and hash tables, so you also need:
 *   #define STB_DS_IMPLEMENTATION 
 * Probably in the future the good idea will be to include parts of stb_ds into this lib
 * directly or put here more basic implementation of hash tables and implement dynamic
 * arrays using realloc.
 *
 * QUICK NOTES:
 *   
 *   Should be used for fast rectangle-rectangle collision detection, for instance in
 *   broadphase scenario or in more simple cases directly.
 * 
 * DOCUMENTATION
 *
 * Limitations:
 *   Key used in a spatial hash is plain int, because we store 2 coordinates in the same
 *   key the size limit is bounded to half size of int, i.e. if int is 32 bits, 
 *   coordinates (after dividing by cell size) are limited to [-65535, 65535] (16 bits).
 *
 *
 * Initialization:
 *   struct isls2d sh;
 *   isls2d_init(&sh, 32.0f, 32.0f);
 *
 * Deinit:
 *   isls2d_clear(&sh);
 *
 * Insert new entity and acquire the id for it:
 *   int id = isls2d_insert(&sh, x, y, width, height, userdata);
 *
 * Remove entity:
 *   isls2d_remove(&sh, id);
 *
 * Update entity:
 *   isls2d_update(&sh, id, new_x, new_y, new_width, new_height);  
 *
 *
 * Additional compilation defines:
 *   ISL_SPATIAL2D_STATIC - static compilation
 *   ISL_SPATIAL2D_DOUBLE - use doubles instead of floats
 *
 * LICENSE
 *
 *   See end of file for license information.
 *
 */

#include "stb_ds.h"

#ifndef ISLS2D_DEF
#ifdef ISL_SPATIAL2D_STATIC
#define ISLS2D_DEF static
#else
#define ISLS2D_DEF extern
#endif
#endif

#define ISLS2D_XMULT (1 << (4*sizeof(int)))
#define ISLS2D_KEY(x, y) ((x)*ISLS2D_XMULT + (y))
#define ISLS2D_X(key) ((x)/ISLS2D_XMULT)
#define ISLS2D_Y(key) ((y)%ISLS2D_XMULT)

#ifndef ISL_SPATIAL2D_DOUBLE
#define isls2d_float  float
#define isls2d__floor floorf
#define isls2d__ceil  ceilf
#else
#define isls2d_float  double
#define isls2d__floor floor
#define isls2d__ceil  ceil
#endif

struct isls2d_entity {
	int id;
	isls2d_float x;
	isls2d_float y;
	isls2d_float width;
	isls2d_float height;
	int xmin;
	int xmax;
	int ymin;
	int ymax;
	const void *data;
	struct {int key; bool value;} *overlaps;
};

struct isls2d {
	struct {int key; int *value;} *cells;
	struct isls2d_entity *entities;
	int *reusable_ids;
	isls2d_float inv_cell_width;
	isls2d_float inv_cell_height;
	bool track_overlap;
};

#ifdef __cplusplus
extern "C" {
#endif

ISLS2D_DEF void isls2d_init(struct isls2d *sh, isls2d_float cell_width, isls2d_float cell_height);
ISLS2D_DEF void isls2d_clear(struct isls2d *sh);
ISLS2D_DEF int isls2d_insert(struct isls2d *sh, isls2d_float x, isls2d_float y, isls2d_float width, isls2d_float height, const void *data);
ISLS2D_DEF void isls2d_remove(struct isls2d *sh, int id);
ISLS2D_DEF void isls2d_update(struct isls2d *sh, int id, isls2d_float x, isls2d_float y, isls2d_float width, isls2d_float height);
#define isls2d_overlaps(x1,y1,w1,h1,x2,y2,w2,h2) ((x1)+(w1)>(x2)&&(x2)+(w2)>(x1)&&(y1)+(h1)>(y2)&&(y2)+(h2)>(y1))

#ifdef __cplusplus
}
#endif
#endif // ISL_INCLUDE_SPATIAL2D_H_

#ifdef ISL_SPATIAL2D_IMPLEMENTATION
#ifndef ISL_SPATIAL2D_IMPLEMENTATION_ONCE 
#define ISL_SPATIAL2D_IMPLEMENTATION_ONCE
#else
#error "ISL_SPATIAL2D_IMPLEMENTATION should be defined once"
#endif // ISL_SPATIAL2D_IMPLEMENTATION_ONCE

#include <math.h>

static int *isls2d__arrsorted_put_if_absent(int *a, int v);
static int *isls2d__arrsorted_del(int *a, int v);
static void isls2d__insert_entity_into_cells(struct isls2d *sh, struct isls2d_entity *e);
static void isls2d__remove_entity_from_cells(struct isls2d *sh, struct isls2d_entity *e);


int *isls2d__arrsorted_put_if_absent(int *a, int v) {
	int n = arrlen(a), l = 0, r = n - 1, m = (l + r) >> 1;
	while (l <= r) {
		if (a[m] < v)      l = m + 1;
		else if (a[m] > v) r = m - 1;
		else return a;
		m = (l + r) >> 1;
	}
	arrins(a, l, v);
	return a;
}

int *isls2d__arrsorted_del(int *a, int v) {
	int n = arrlen(a), l = 0, r = n - 1, m = (l + r) >> 1;
	while (l <= r) {
		if (a[m] < v)      l = m + 1;
		else if (a[m] > v) r = m - 1;
		else {
			arrdel(a, m);
			return a;
		}
		m = (l + r) >> 1;
	}
	return a;
}

void isls2d__insert_entity_into_cells(struct isls2d *sh, struct isls2d_entity *e) {
	int id = e->id, xmin = e->xmin, xmax = e->xmax, ymin = e->ymin, ymax = e->ymax;
	for (int x = e->xmin; x < xmax; x++) {
		for (int y = ymin; y < ymax; y++) {
			int key = ISLS2D_KEY(x, y);
			int *cell_ids = hmget(sh->cells, key);
			int n = arrlen(cell_ids);
			if (track_overlap) {
				for (int i = 0; i < n; i++) {
					struct isls2d_entity *o = &sh->entities[i];
					if (isls2d_overlaps(e->x, e->y, e->width, e->height, o->x, o->y, o->width, o->height)) {
						isls2d__arrsorted_put_if_absent(e->overlaps, o->id);
						isls2d__arrsorted_put_if_absent(o->overlaps, e->id);
						//hmput(e->overlaps, o->id, true);
						//hmput(o->overlaps, e->id, true);
					}
				}
			}
			arrput(cell_ids, id);
			if (n == 1) hmput(sh->cells, key, cell_ids);
		}
	}
}

void isls2d__remove_entity_from_cells(struct isls2d *sh, int id, int xmin, int xmax, int ymin, int ymax) {
	for (int x = xmin; x < xmax; x++) {
		for (int y = ymin; y < ymax; y++) {
			int key = ISLS2D_KEY(x, y);
			int *cell_ids = hmget(sh->cells, key);
			int n = arrlen(cell_ids);
			for (int i = 0; i < n; i++) {
				if (ids[i] == id) {
					arrdelswap(cell_ids, i);
					if (track_overlap) {
						for (int j = 0; j < n-1; j++) {
							struct isls2d_entity *o = &sh->entities[i];
							isls2d__arrsorted_del(e->overlaps, o->id);
							isls2d__arrsorted_del(o->overlaps, e->id);
							//hmdel(e->overlaps, o->id);
							//hmdel(o->overlaps, e->id);
						}
					}
					break;		
				}
			}
			if (arrlen(cell_ids) == 0) {
				arrfree(ids);
				(void)hmdel(sh->cells, key);
			}
		}
	}
}

void isls2d_init(struct isls2d *sh, isls2d_float cell_width, isls2d_float cell_height) {
	*sh = (struct isls2d) {NULL, NULL, NULL, 1.0f / cell_width, 1.0f / cell_height};
}

void isls2d_clear(struct isls2d *sh) {
	int n = hmlen(sh->cells);
	for (int i = 0; i < n; i++) {
		arrfree(sh->cells[i].value);
	}
	hmfree(sh->cells);
	n = arrlen(sh->entities);
	for (int i = 0; i < n; i++) {
		//hmfree(sh->entities[i].overlaps);
		arrfree(sh->entities[i].overlaps);
		sh->entities[i].overlaps = NULL;
	}
	arrfree(sh->entities);
	arrfree(sh->reusable_ids);
	sh->cells = NULL;
	sh->entities = NULL;
	sh->reusable_ids = NULL;
}

int isls2d_insert(struct isls2d *sh, isls2d_float x, isls2d_float y, isls2d_float width, isls2d_float height, const void *data) {
	int id;
	int xmin = isls2d__floor(x * sh->inv_cell_width);
	int xmax = isls2d__ceil((x + width) * sh->inv_cell_width);
	int ymin = isls2d__floor(y * sh->inv_cell_height);
	int ymax = isls2d__ceil((y + height) * sh->inv_cell_height);
	struct isls2d_entity entity = (struct isls2d_entity) {-1, x, y, width, height, xmin, xmax, ymin, ymax, data, NULL};
	if (arrlen(sh->reusable_ids) > 0) {
		entity.id = arrpop(sh->reusable_ids);
		sh->entities[entity.id] = entity;
	} else {
		entity.id = arrlen(sh->entities);
		arrpush(sh->entities, entity);
	}
	isls2d__insert_entity_into_cells(sh, id, xmin, xmax, ymin, ymax);
	return id;
}

void isls2d_remove(struct isls2d *sh, int id) {
	if (id >= arrlen(sh->entities)) return;
	struct isls2d_entity entity = sh->entities[id];
	if (entity.id != id) return;
	isls2d__remove_entity_from_cells(sh, id, entity.xmin, entity.xmax, entity.ymin, entity.ymax);
	sh->entities[id].id = -1;
	int n = arrlen(sh->entities[id].overlaps);
	for (int i = 0; i < n; i++) {
		isls2d__arrsorted_del(sh->entities[i].overlaps, id);
	}
	arrfree(sh->entities[id].overlaps);
	sh->entities[id].overlaps = NULL;
	if (id == arrlen(sh->entities)) {
		(void)arrpop(sh->entities);
	} else {
		(void)arrput(sh->reusable_ids, id);
	}
}

void isls2d_update(struct isls2d *sh, int id, isls2d_float x, isls2d_float y, isls2d_float width, isls2d_float height) {
	if (id >= arrlen(sh->entities)) return;
	struct isls2d_entity entity = sh->entities[id];
	if (entity.id != id) return;
	int xmin = isls2d__floor(x * sh->inv_cell_width);
	int xmax = isls2d__ceil((x + width) * sh->inv_cell_width);
	int ymin = isls2d__floor(y * sh->inv_cell_height);
	int ymax = isls2d__ceil((y + height) * sh->inv_cell_height);
	if (sh->track_overlap || xmin != entity.xmin || xmax != entity.xmax || ymin != entity.ymin || ymax != entity.ymax) {
		isls2d__remove_entity_from_cells(sh, id, entity.xmin, entity.xmax, entity.ymin, entity.ymax);
		isls2d__insert_entity_into_cells(sh, id, xmin, xmax, ymin, ymax);
	}
	sh->entities[id] = (struct isls2d_entity) {entity.id, x, y, width, height, xmin, xmax, ymin, ymax, entity.data, entity.overlaps};
}

/*
------------------------------------------------------------------------------
This software is available under 2 licenses -- choose whichever you prefer.
------------------------------------------------------------------------------
ALTERNATIVE A - MIT License
Copyright (c) 2022 Ilya Kolbin
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------
ALTERNATIVE B - Public Domain (www.unlicense.org)
This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
software, either in source code form or as a compiled binary, for any purpose,
commercial or non-commercial, and by any means.
In jurisdictions that recognize copyright laws, the author or authors of this
software dedicate any and all copyright interest in the software to the public
domain. We make this dedication for the benefit of the public at large and to
the detriment of our heirs and successors. We intend this dedication to be an
overt act of relinquishment in perpetuity of all present and future rights to
this software under copyright law.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
------------------------------------------------------------------------------
*/
#endif // ISL_SPATIAL2D_IMPLEMENTATION
