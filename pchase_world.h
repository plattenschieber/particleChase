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
        double              t;  /* current time */
        double              delta_t;    /* timestep length */
        double              t_end;      /* end time */
        unsigned int        n_particles;        /* number of particles */
        unsigned int        step;       /* current step */
        double              length[DIM];        /* length */
        p4est_t            *p4est;      /* a pointer to the allocated p4est */
}
                    pchase_world_t;

/** initialize a world with static parameter
 *
 * \param[in] p4est a forest we want to integrate into our world
 * \return	        a ready to use world
 */
pchase_world_t     *pchase_world_init(p4est_t * p4est);

/** start the simulation
 *
 * \param [in] W	the world we are working on
 */
void                pchase_world_simulate(pchase_world_t * W);

/** update positions
 *
 * \param [in] W	the world we are working on
 */
void                pchase_world_update_x(pchase_world_t * W);

/** generate a particle via random distribution inside the pchase_worlds boundaries
 *
 * \param [in] W	the world we are working on
 * \return          generated particle inside the world
 */
pchase_particle_t  *pchase_world_random_particle(pchase_world_t * W);

/** insert a particle into its belonging quadrant
 *
 * \param [in] W	the world into which we are operating
 * \param [in] p	the particle to be inserted
 * \return          return the mini quad where the particle stays
 */
p4est_quadrant_t   *pchase_world_insert_particle(pchase_world_t * W, pchase_particle_t * p);

/** prints out all particles into a XYZ file
 *
 * \param [in] W	the world which shall be printed
 */
void                pchase_world_print_particlesXYZ(pchase_world_t * W);

/** take out a particle of its world and find out where it stays in p4est
 *
 * \param [in] W	the world we are working on
 * \param [in] p	a particle lying inside the world
 * \return 		    the quadrant to which the given particle belongs
 */
p4est_quadrant_t   *pchase_translate_particle_to_p4est(pchase_world_t * W, pchase_particle_t * p);

/** find the owner of a given quadrant
 *
 * \param [in] q	quadrant whos owner rank shall be found
 * \return 		    owner rank
 */
int                 pchase_quadrant_is_in_proc(p4est_quadrant_t * q);

/** callback indicating if a given mini quad (*point), which encloses a
 * particle, belongs to the current quad set marker to this quad
 * if (current quad is enclosing particle && this.quad.level > marker.level
 *
 * \param[in] p4est     the forest we are working on
 * \param[in] which_tree the tree we are looking at
 * \param[in] quadrant  current quadrant
 * \param[in] is_leaf   is this a leaf or not
 * \param[in] point     our mini quad, enclosing the inserted particle
 * \return              true if we found a quadrant that is enclosing our mini quad
 */
static int
search_fn(p4est_t * p4est, p4est_topidx_t which_tree,
          p4est_quadrant_t * quadrant, int is_leaf, void *point);
#endif
