#include "pchase_world.h"

pchase_world_t *pchase_world_init(FILE *parameter, FILE *particles) {
	int i;
	/* get some place for the pchase_world */
	pchase_world_t *W;
	if (!(W=malloc(sizeof(pchase_world_t)))) {
		fprintf(stderr, "Could not allocate memory for the World\n");
		exit(EXIT_FAILURE); 
	}
	/* set all parameters */
	W->t 			= 0.0;
	W->delta_t 		= 0.1;
	W->t_end		= 1.0;
	W->n_particles 	= 0;
	W->step 		= 0;
	for (i=0; i<DIM; i++)
		W->length[i] = 2.0;
	/* reset seed */
	srand(time(NULL));
	return W;
}

void pchase_world_simulate(pchase_world_t *W) {
	/* simulate until the end has come */
	while (W->t <= W->t_end) {
		pchase_world_update_x(W);
#ifdef DEBUG
		printf("Actual time on earth: %f\n",W->t);
#endif
		W->t+=W->delta_t;
		W->step++;
	}
}

void pchase_world_update_x(pchase_world_t *W) {
	/* evaluate potential */
	printf("EVALUATE VELOCITY FIELD - NOT IMPLEMENTED YET");
}


pchase_particle_t * pchase_world_random_particle(pchase_world_t *W) {
	int i;
	pchase_particle_t *p = malloc(sizeof(pchase_particle_t));
	for (i=0; i<DIM; i++)
		p->x[i] = W->length[i] * rand()/(RAND_MAX + 1.);
#ifdef DEBUG
	p->ID = W->n_particles;
#endif 
	W->n_particles++; 

	return p;
}

void pchase_world_print_particlesXYZ(pchase_world_t *W){
	prinft("PRINT XYZ - NOT IMPLEMENTED YET");
}
