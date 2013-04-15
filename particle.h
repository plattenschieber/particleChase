#ifndef PARTICLE_H
#define PARTICLE_H
#include <stdlib.h>
#include <stdio.h>
#define DIM 2


typedef struct {
	int ID;
	double x[DIM];
	double v[DIM];
} PARTICLE;

typedef struct NODE {
	void *data;
	struct NODE *next;
} NODE;

typedef struct {
	NODE *first;
} LIST;

void list_free(LIST * L);
/* delete the 'cursor' from the list and leave everything clean */
void list_delete(LIST *L, NODE *cursor);

#endif
