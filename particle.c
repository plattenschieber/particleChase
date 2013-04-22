#include "particle.h"
#include "world.h"
#include <stdlib.h>
#include <stdio.h>

NODE *list_create(LIST *L, PARTICLE *particle) {
	NODE *node;
	if( !(node=malloc(sizeof(NODE))) ) return NULL;
	node->particle = particle;
	node->next = NULL;
	return node;
}

NODE *list_insert(LIST *L, NODE *cursor, PARTICLE *particle) {
	NODE *newnode;
	/* create 'newnode' only if there is enough space */ 
	if( (newnode=list_create(L, particle)) ) {
		/* when 'cursor' appears in the list */
		if (cursor != L->last) {
			newnode->next = cursor->next;
			cursor->next = newnode;
		}
		/* 'cursor' is either NULL or 'L->last', so append it to the end */
		else {
			newnode->next = NULL;
			/* newnode is firstly new 'next' of 'L->last' and than it's 'L->last' (we're reading from right to left) */
			L->last = L->last->next = newnode;
		return newnode;
	}
	/* we can't continue */
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
        free(it->particle);
        free(it);
    }
    free(L);
}

