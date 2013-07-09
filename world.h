#ifndef WORLD_H
#define WORLD_H
#include "particle.h"
#include <stdio.h> 
typedef struct {
	/* a list of cells with particle lists */
	LIST *cells;
	// current time */
	double t;
	// timestep */
	double delta_t;
	// end time */
	double t_end;
	// number of overall particles */
	unsigned int n_particles;
	/* actual step */
	unsigned int step;
} WORLD;

/* start the simulation */
void simulate(LIST *L);
/* update the positions */
void update_x(LIST *L);
/* initialize the world's parameters */
void world_init(FILE *F);

#endif
