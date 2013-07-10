#ifndef WORLD_H
#define WORLD_H
#include "list.h"
#include <stdio.h> 
typedef struct {
	/* a list of cells with particle lists */
	LIST *cells;
	// current time */
	double t;
	// timestep length */
	double delta_t;
	// end time */
	double t_end;
	// number of overall particles */
	unsigned int n_particles;
	/* current step */
	unsigned int step;
	/* length */ 
	double length[DIM]; 
} WORLD;

/* initialize the world's parameters */
WORLD *world_init(FILE *F);
/* start the simulation */
void simulate(WORLD *W);
/* update positions */
void update_x(WORLD *W);
/* update velocity */
void update_v(WORLD *W);

#endif
