#ifndef PARTICLE_H
#define PARTICLE_H
#define DIM 2

typedef struct {
	int ID;
	double x[DIM];
	double v[DIM];
} PARTICLE;

typedef struct NODE {
	PARTICLE *particle;
	struct NODE *next;
} NODE;

typedef struct {
	NODE *first;
	NODE *last;
} LIST;

/* get some space for a new 'LIST' and initialize it with NULL */
LIST *list_init();
/* create a new node with a 'particle' and return its adress - only intern use */
NODE *list_create_node(LIST *L, PARTICLE *particle);
/* insert a new new node with a 'particle' behind 'cursor' and set all pointers accordingly */
NODE *list_insert(LIST *L, NODE *cursor, PARTICLE *particle);
/* delete the 'cursor' from the list and leave everything clean */
void list_delete(LIST *L, NODE *cursor);
/* delete first element of the list and leave everything clean */
void list_pop_first(LIST *L);
/* delete every element of the list and free it */
void list_free(LIST *L);
/* print every particle */
void print_particles(LIST *L);
#endif
