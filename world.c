#include "world.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

WORLD *world_init(FILE *parameter, FILE *particles) {
	WORLD *W;
	if (!(W=malloc(sizeof(WORLD)))) {
	/* get some place for the pchase_world */
		fprintf(stderr, "Could not allocate memory for the World\n");
		exit(EXIT_FAILURE); 
	}
	/* set all parameters */
	W->n_particles = 0;
	W->step = 0;
	world_read_parameter(W, parameter);
	/* world_read_particles(W, particles); */
	/* reset seed */
	srand(time(NULL));
	world_randomfill(W);
	/* fill in some particles randomly */
	return W;
}

void world_simulate(WORLD *W) {
	/* simulate until the end has come */
	while (W->t <= W->t_end) {
		world_update_x(W);
#ifdef DEBUG
		printf("Actual time on earth: %f\n",W->t);
#endif
		print_particles(W->particles);
		W->t+=W->delta_t;
		W->step++;
	}
}

}


PARTICLE * random_particle() {
	int i;
	PARTICLE *p = malloc(sizeof(PARTICLE));
	for (i=0; i<DIM; i++)
		p->x[i] = rand()/(RAND_MAX + 1.);
	p->ID = PARTICLE_COUNT;
	PARTICLE_COUNT++; 
	return p;
}
PARTICLE * world_random_particle(WORLD *W) {
	int i;
	PARTICLE *p = malloc(sizeof(PARTICLE));
	for (i=0; i<DIM; i++)
		p->x[i] = W->length[i] * rand()/(RAND_MAX + 1.);
	p->ID = W->n_particles;
	list_insert(W->particles, NULL, p);
	W->n_particles++; 
	return p;
}

}
