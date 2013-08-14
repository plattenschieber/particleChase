#ifndef WORLD_H
#define WORLD_H
#define DEBUG
#define DIM 2

#include <stdio.h> 
/* pchase_world_t holds the entire information of our simulation */
typedef struct {
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
WORLD *world_init(FILE *parameter, FILE *particles);
/* start the simulation */
void world_simulate(WORLD *W);
/* update positions */
void world_update_x(WORLD *W);
/* return a randomly assigned particle inside the world */
PARTICLE * random_particle();
/* return a particle via random distribution inside the pchase_worlds boundaries */
/* print XYZ file */
void world_print_particlesXYZ(WORLD *W);
#endif
