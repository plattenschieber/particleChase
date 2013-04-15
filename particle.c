#include "particle.h"

NODE *list_create(LIST *L, void *data) {
	NODE *node;
	if( !(node=malloc(sizeof(NODE))) ) return NULL;
	node->data=data;
	node->next=NULL;
	return node;
}

NODE *list_insert(LIST *L, NODE *cursor, void *data) {
	NODE *newnode;
	if( (newnode=list_create(L, data)) ) {
        newnode->next = cursor->next;
        cursor->next = newnode;
		return newnode;
	}
	else {
		fprintf(stderr, "Could not allocate memory for new node\n");
        exit(EXIT_FAILURE); 
	}
}

void list_delete(LIST *L, NODE *cursor) {
	NODE *it = L->first;
	while(it->next && it->next!=cursor) it=it->next;
	/* ToDo */
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

