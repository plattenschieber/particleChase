#ifndef PCHASE_WORLD_H
#define PCHASE_WORLD_H
#define DEBUG
#define DIM 2

#include "pchase_particle.h"
#include <stdlib.h>

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
} pchase_world_t;

/* initialize the pchase_world's parameters */
pchase_world_t *pchase_world_init();
/* start the simulation */
void pchase_world_simulate(pchase_world_t *W);
/* update positions */
void pchase_world_update_x(pchase_world_t *W);
/* return a particle via random distribution inside the pchase_worlds boundaries */
pchase_particle_t * pchase_world_random_particle(pchase_world_t *W);
/* print XYZ file */
void pchase_world_print_particlesXYZ(pchase_world_t *W);
#endif
