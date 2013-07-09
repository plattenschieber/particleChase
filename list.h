#ifndef LIST_H
#define LIST_H
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

/* get some space for a new 'LIST' and initialize it with NULL */
LIST *list_init();
/* create a new node with 'data' and return its adress - only intern use */
NODE *list_create_node(LIST *L, void *data);
/* insert a new new node with 'data' behind 'cursor' and set all pointers accordingly */
NODE *list_insert(LIST *L, NODE *cursor, void *data);
/* delete the 'cursor' from the list and leave everything clean */
void list_delete(LIST *L, NODE *cursor);
/* delete first element of the list and leave everything clean */
void list_pop_first(LIST *L);
/* delete every element of the list and free it */
void list_free(LIST *L);
/* print every particle */
void print_particles(LIST *L);
#endif
