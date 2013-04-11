#include <stdio.h>
#include <stdlib.h>
#include "particle.h"

int main() {
	int i;
	PARTICLE particle;
	NODE *p;
	LIST *L;
	particle.x[0] = 1.3; particle.x[1] = 2.3;
	particle.ID = 2;

	/* loop through the list L */
	//for (p = L->first; p; p = p->next) {
		/* code */
	

	for (i = 0; i < DIM; i++) {
		printf("%f\n",particle.x[i]);
	}
	printf("%i\n",particle.ID);
	return 0;
}
