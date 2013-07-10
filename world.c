#include "world.h"
#include <stdio.h>
#include <stdlib.h>

WORLD *world_init(FILE *parameter, FILE *particles) {
	/* get some place for the world */
	WORLD *W;
	if (!(W=malloc(sizeof(WORLD)))) {
		fprintf(stderr, "Could not allocate memory for the World\n");
		exit(EXIT_FAILURE); 
	}
	/* get some place for its cells */
	W->cells = list_init();
	/* read in and set all parameters */
	W->n_particles = 0;
	W->step = 0;
	world_read_particles(W, particles);
	world_read_parameter(W, parameter);
	return W;
}

void world_simulate(WORLD *W) {
	/* simulate until the end has come */
	while (W->t <= W->t_end) {
		world_update_x(W);
		printf("Actual time on earth: %f\n",W->t);
		W->t+=W->delta_t;
		W->step++;
	}
}

void world_update_x(WORLD *W) {
	NODE *it;
	PARTICLE *p;
	unsigned int i;
 	for (it = W->cells->first; it != NULL; it = it->next)
		for (i=0; i<DIM; i++) {
			/* cast data to particle */
			p = it->data;
			/* do some calculation (move more in y coordinate than x)*/
			p->x[i] += 0.112+i/2;
		}
}

void world_update_v(WORLD *W) {
	return;
}

void world_read_particles(WORLD *W, FILE *particles) {
	int ID;
	double x[DIM];
	/* read in particles in this form: ID x[0] x[1] */
	while (fscanf(particles,"%i %lf %lf",&ID, &x[0], &x[1]) == 3) { 
		PARTICLE *tmp = malloc(sizeof(PARTICLE));
		tmp->ID = ID;
		printf("%i\t",tmp->ID);
		int i;
		for (i = 0; i < DIM; i++){
			tmp->x[i] = x[i];
			printf("%f\t",tmp->x[i]);
		}
		printf("\n");

		list_insert(W->cells, NULL, tmp);
		W->n_particles++; 
		printf("particle inserted\n");
	}
	printf("Read %i Particle(s)\n", W->n_particles);
	
}
void world_read_parameter(WORLD *W, FILE *parameter) {
	char *option;
	int tmp1;
	double tmp2;
	printf("read parameter\n");
	fscanf(parameter,"%s", option);
	printf("%s",option);
}
