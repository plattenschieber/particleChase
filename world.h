#ifndef WORLD_H
#define WORLD_H
#define DEBUG
#define DIM 2

#include <stdio.h> 
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
/* update velocity */
void world_update_v(WORLD *W);

/* read in particles and save them into their corresponding cell */
void world_read_particles(WORLD *W, FILE *F);
/* read world parameter from file */
void world_read_parameter(WORLD *W, FILE *F);
/* fill in MAX_PARTICLES particles */
void world_randomfill(WORLD *W);
/* fill the world with random particles */
PARTICLE * world_random_particle(WORLD *W);
/* return a randomly assigned particle inside the world */
PARTICLE * random_particle();
/* print XYZ file */
void world_print_particlesXYZ(WORLD *W);
#endif
