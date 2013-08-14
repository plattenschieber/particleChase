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
	double 		t;		/* current time */
	double 		delta_t;	/* timestep length */
	double 		t_end;		/* end time */
	unsigned int 	n_particles;	/* number of overall particles */
	unsigned int 	step;		/* current step */
	double 		length[DIM];	/* length */ 
	p4est_t		p4est;		/* a pointer to the allocated p4est */
} pchase_world_t;

/** initialize a world with static parameter
 * \return		a ready to use world
 */
pchase_world_t *pchase_world_init(p4est_t p4est);

/** start the simulation 
 * \param [in] W	the world we are working on
 */
void pchase_world_simulate(pchase_world_t *W);

/** update positions
 * \param [in] W	the world we are working on
 */
void pchase_world_update_x(pchase_world_t *W);

/** return a particle via random distribution inside the pchase_worlds boundaries
 * \param [in] W	the world we are working on
 */
pchase_particle_t * pchase_world_random_particle(pchase_world_t *W);

/** insert a particle into its belonging quadrant and return the quadrant
 * \param [in] W	the world into which we are operating
 * \param [in] p	the particle to be inserted 
 */
p4est_quadrant_t * pchase_world_insert_particle(pchase_world_t *W, pchase_particle_t *p);

/** prints out all particles into a XYZ file
 * \param [in] W	the world which shall be printed 
 */
void pchase_world_print_particlesXYZ(pchase_world_t *W);

/** converts a worlds (particle) coordinate into a p4est quadrant coordinate 
 * \param [in] coord	world coordinate to be converted
 */
p4est_quadrant_t * pchase_translate_particle_to_p4est(pchase_particle_t *p);
#endif
