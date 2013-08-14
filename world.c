#include "pchase_world.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

pchase_pchase_world_t *pchase_world_init(FILE *parameter, FILE *particles) {
	/* get some place for the pchase_world */
	pchase_pchase_world_t *W;
	if (!(W=malloc(sizeof(pchase_pchase_world_t)))) {
		fprintf(stderr, "Could not allocate memory for the World\n");
		exit(EXIT_FAILURE); 
	}
	/* set all parameters */
	W->n_particles = 0;
	W->step = 0;
	/* reset seed */
	srand(time(NULL));
	/* fill in some particles randomly */
	pchase_world_randomfill(W);
	return W;
}

void pchase_world_simulate(pchase_pchase_world_t *W) {
	/* simulate until the end has come */
	while (W->t <= W->t_end) {
		pchase_world_update_x(W);
#ifdef DEBUG
		printf("Actual time on earth: %f\n",W->t);
#endif
		print_particles(W->particles);
		W->t+=W->delta_t;
		W->step++;
	}
}

void pchase_world_update_x(pchase_pchase_world_t *W) {
	/* evaluate potential */
	prinft("EVALUATE VELOCITY FIELD - NOT IMPLEMENTED YET");
}


pchase_particle_t * pchase_world_random_particle(pchase_pchase_world_t *W) {
	int i;
	pchase_particle_t *p = malloc(sizeof(pchase_particle_t));
	for (i=0; i<DIM; i++)
		p->x[i] = W->length[i] * rand()/(RAND_MAX + 1.);
#ifdef DEBUG
	p->ID = pchase_world_pcounter;
#endif 
	pchase_world_pcounter++; 

	return p;
}

void pchase_world_print_particlesXYZ(pchase_world_t *W){
	prinft("PRINT XYZ - NOT IMPLEMENTED YET");
}
