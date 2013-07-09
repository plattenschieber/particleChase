#include "world.h"
#include "world.h"
#include <stdio.h>
#include <stdlib.h>

	}
}
void world_init(FILE *F) {
	/* read from File instead */
	world_t = 0.0;
	world_t_end = 2.0;
	world_delta_t = 0.1;
void simulate(WORLD *W) {
	while(W->t < W->t_end) {
		update_x(W);
		printf("Actual time on earth: %f\n",W->t);
		W->t+=W->delta_t;
		W->step++;
	}
}

void update_x(WORLD *W) {
	/* do some stuff */
	NODE *it;
	unsigned int i;
 	for (it = W->cells->first; it != NULL; it = it->next)
		for(i=0; i<DIM; i++)
			/* do some calculation (move more in y coordinate than x)*/
			it->particle->x[i] += 0.112+i/2;
}

