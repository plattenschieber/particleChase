#ifndef PARTICLE_H
#define PARTICLE_H
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
	NODE *last;
} LIST;

/* create a new node with 'data' and return its adress */
NODE *list_create(LIST *L, void *data);
/* insert a new node with 'data' behind 'cursor' and reset 'L->first' if 'cursor==NULL' */
NODE *list_insert(LIST *L, NODE *cursor, void *data);
/* delete the 'cursor' from the list and leave everything clean */
void list_delete(LIST *L, NODE *cursor);
/* remove the first element of the list and leave everything clean */
void list_pop_first(LIST *L);
/* delete every element of the list and free it */
void list_free(LIST *L);

#endif
