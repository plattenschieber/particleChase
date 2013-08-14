#ifndef PCHASE_WORLD_H
#define PCHASE_WORLD_H
#define DEBUG
#define DIM 2

#include "pchase_particle.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "p4est_extended.h"

/* pchase_world_t holds the entire information of our simulation */
typedef struct {
	double t;			/* current time */
	double delta_t;			/* timestep length */
	double t_end;			/* end time */
	unsigned int n_particles;	/* number of overall particles */
	unsigned int step;		/* current step */
	double length[DIM];		/* length */ 
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
