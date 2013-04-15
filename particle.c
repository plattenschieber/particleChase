#include "particle.h"

/* insert data right behind cursor. If cursor equals NULL, insert at the beginning. Return inserted node or NULL if an errror occurs */
NODE * list_insert(LIST * L, NODE * cursor, double data) {
	NODE * insert;
	if ( (insert=malloc(sizeof(NODE))) != NULL ) {
		insert->data = &data;
		insert->next = cursor->next;
		cursor->next = insert;
NODE *list_create(LIST *L, void *data) {
	NODE *node;
	if( !(node=malloc(sizeof(NODE))) ) return NULL;
	node->data=data;
	node->next=NULL;
	return node;
}

	}
	return insert;
}
