LIST * list_create() {
	LIST *L;
    if ( (L=malloc(sizeof(LIST))) != NULL ) 
        L->first = NULL;
    return L;
}

/* insert data right behind cursor. If cursor equals NULL, insert at the beginning. Return inserted node or NULL if an errror occurs */
NODE * list_insert(LIST * L, NODE * cursor, double data) {
	NODE * insert;
	if ( (insert=malloc(sizeof(NODE))) != NULL ) {
		insert->data = &data;
		insert->next = cursor->next;
		cursor->next = insert;

	}
	return insert;
}
