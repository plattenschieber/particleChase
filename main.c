#include <stdio.h>
#include <stdlib.h>
#include "particle.h"
#include "world.h"

int main() {
	int i;
	PARTICLE p1;
   	PARTICLE p2;
	LIST *L;
	NODE *it;

	p1.x[0] = 1.3; p1.x[1] = 2.3;
	p1.ID = 2;
	L->first = list_create(L, &p1);

	p2.x[0] = 1.5; p2.x[1] = 4.8;
	p2.ID = 1;
	list_insert(L, L->first, &p2);

	/* loop through the list L */
	for (it = L->first; it; it = it->next) {
		printf("%i\t",it->particle->ID);
		for (i = 0; i < DIM; i++) 
			printf("%f\t",it->particle->x[i]);
		printf("\n");
	}
	simulate(L);
	return 0;
}
