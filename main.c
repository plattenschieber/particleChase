#include <stdio.h>
#include <stdlib.h>
#include "particle.h"
#include "world.h"

int main() {
	int i;
	PARTICLE particle;
	LIST *L;
	NODE *it;

	particle.x[0] = 1.3; particle.x[1] = 2.3;
	particle.ID = 2;
	L->first = list_create(L, &particle);

	particle.x[0] = 1.5; particle.x[1] = 4.8;
	particle.ID = 1;
	list_insert(L, L->first, &particle);

	/* loop through the list L */
	for (it = L->first; it; it = it->next) {
		PARTICLE tmp;
	   	//tmp = &(it->data);
		printf("%i\t",tmp.ID);
		for (i = 0; i < DIM; i++) 
			printf("%f\t",tmp.x[i]);
	}
	simulate(L);
	return 0;
}
