#include "list.h"
#include "world.h"
#include <stdlib.h>
#include <stdio.h>

LIST *list_init() {
	LIST *L;
	if( !(L=malloc(sizeof(LIST))) ) {
		fprintf(stderr, "Could not allocate memory for the List\n");
		exit(EXIT_FAILURE); 
	}
	L->first = NULL;
	L->last = NULL;
	return L; 
}

NODE *list_create_node(LIST *L, void *data) {
	NODE *node;
	if( !(node=malloc(sizeof(NODE))) ) return NULL;
	node->data = data;
	node->next = NULL;
	return node;
}

NODE *list_insert(LIST *L, NODE *cursor, void *data) {
	NODE *newnode;
	/* create 'newnode' only if there is enough space */ 
	if( (newnode=list_create_node(L, data)) ) {
		/* 'cursor' appears in the list */
		if (cursor) {
			newnode->next = cursor->next;
			cursor->next = newnode;
			if (cursor == L->last) L->last = newnode;
		}
		/* 'cursor' is not specified (== NULL) */
		else {
			/* 'L' is nonempty -> append 'newnode' */
			if (L->last) L->last = L->last->next = newnode;
			/* 'L' is empty */
			else L->first = L->last = newnode;
		}
		return newnode;
	}
	/* there is not enough space -> we can't continue */
	else {
		fprintf(stderr, "Could not allocate memory for new node\n");
		exit(EXIT_FAILURE); 
	}
}

void list_delete(LIST *L, NODE *cursor) {
	NODE *it = L->first, *tmp;
	/* error handling */
	if (cursor == L->first) 
		list_pop_first(L);
	/* actual delete function */
	else {
		/* find the cursor in the list */
		while(it->next && it->next!=cursor) it=it->next;

		/* if cursor was 'L->last', reset it */
		if (cursor == L->last) 
			L->last = it;

		/* it wasn't 'L->last' so we need to free something */
		else {
			/* save 'cursor->next' and free 'cursor' */
			tmp = it->next->next;
			free(it->next->data);
			free(it->next);
			/* reset pointer accordingly */	
			it->next = tmp;
		}
	}
}

void list_pop_first(LIST *L) {
	/* save second item and delete the first */
	NODE *tmp = L->first->next;
	free(L->first->data);
	free(L->first);
	L->first = tmp;
} 

void list_free(LIST *L) {
	NODE *it, *tmp;
	for (it = L->first; it != NULL; it = tmp) {
		tmp = it->next;
		free(it->data);
		free(it);
	}
	free(L);
}

void print_particles(LIST *L) {
	int i;
	NODE *it;
	for (it = L->first; it; it = it->next) {
		printf("%i\t",((PARTICLE*)(it->data))->ID);
		for (i = 0; i < DIM; i++) 
			printf("%f\t",((PARTICLE*)(it->data))->x[i]);
		printf("\n");
	}
}

