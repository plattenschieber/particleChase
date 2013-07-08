#include "world.h"
#include "world.h"
#include <stdio.h>
#include <stdlib.h>

void simulate(LIST *L) {
	while(world_t < world_t_end) {
		update_x(L);
		printf("Actual time on earth: %f\n",world_t);
		world_t+=world_delta_t;
		world_step++;
	}
}
void world_init(FILE *F) {
	/* read from File instead */
	world_t = 0.0;
	world_t_end = 2.0;
	world_delta_t = 0.1;
}

void update_x(LIST *L) {
	/* do some stuff */
	NODE *it;
	unsigned int i;
 	for (it = L->first; it != NULL; it = it->next)
		for(i=0; i<DIM; i++)
			/* do some calculation (move more in y coordinate than x)*/
			it->particle->x[i] += 0.112+i/2;
}

