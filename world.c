#include "world.h"
#include <stdio.h>
#include <stdlib.h>

WORLD *world_init(FILE *F) {
	WORLD *W;
	if( !(W=malloc(sizeof(WORLD))) ) {
		fprintf(stderr, "Could not allocate memory for the World\n");
		exit(EXIT_FAILURE); 
	}
	/* read from File instead */
	W->t = 0.0;
	W->t_end = 2.0;
	W->delta_t = 0.1;
	W->cells = list_init();
	W->step = 0;
	W->n_particles = 0;
	return W;
}

void simulate(WORLD *W) {
	while(W->t < W->t_end) {
		update_x(W);
		printf("Actual time on earth: %f\n",W->t);
		W->t+=W->delta_t;
		W->step++;
	}
}

void update_x(WORLD *W) {
	NODE *it;
	PARTICLE *p;
	unsigned int i;
 	for (it = W->cells->first; it != NULL; it = it->next)
		for(i=0; i<DIM; i++) {
			/* cast data to particle */
			p = it->data;
			/* do some calculation (move more in y coordinate than x)*/
			p->x[i] += 0.112+i/2;
		}
}

void update_v(WORLD *W){
	return;
}
