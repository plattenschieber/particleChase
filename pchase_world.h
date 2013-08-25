#ifndef PCHASE_WORLD_H
#define PCHASE_WORLD_H
#define DEBUG
#define PRINTGNUPLOT
/* #define PRINTXYZ */
#define DIM 2

#include "pchase_particle.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "p4est_extended.h"
#include "p4est_search.h"
#include "p4est_iterate.h"
#include "p4est_communication.h"
#include "sc_notify.h"
#include "p4est_vtk.h"

/* pchase_world_t holds the entire information of our simulation */
typedef struct {
        double              t;  /* current time */
        double              delta_t;    /* timestep length */
        double              t_end;      /* end time */
        unsigned int        n_particles;        /* number of overall
                                                 * particles */
        unsigned int        step;       /* current step */
        double              length[DIM];        /* length */
        p4est_t            *p4est;      /* a pointer to the allocated p4est */
        p4est_init_t        init_fn;
        p4est_refine_t      refine_fn;
        p4est_coarsen_t     coarsen_fn;
        p4est_replace_t     replace_fn;
        p4est_search_query_t search_fn;
        p4est_iter_volume_t viter_fn;
        p4est_iter_volume_t destroy_fn;
        p4est_iter_volume_t print_fn;
        p4est_iter_volume_t update_x_fn;
        sc_list_t          *particle_push_list;
        sc_array_t         *particles_to;
}
                    pchase_world_t;

/** initialize a world with static parameter
 *
 * \return	        a ready to use world
 */
pchase_world_t     *pchase_world_init();

/** initialize everything depending on p4est
 *
 * \param[in] the world we are acting on
 * \param[in] p4est a forest we want to integrate into our world
 */
void                pchase_world_init_p4est(pchase_world_t * W, p4est_t * p4est);

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

/** insert all particles int the particle_push_list into their belonging quadrants
 *
 * \param [in] W	the world into which we are operating
 */
void                pchase_world_insert_particles(pchase_world_t * W);

/** prints out all particles into a XYZ file
 *
 * \param [in] W	the world which shall be printed
 */
void                pchase_world_print_particlesXYZ(pchase_world_t * W);

/** take out a particle of its world and find out where it stays in p4est
 *
 * \param [in] W	the world we are working on
 * \param [in] p	a particle lying inside the world
 * \param [out]q    the quadrant to which the given particle belongs
 */
void                pchase_translate_particle_to_p4est(pchase_world_t * W, const pchase_particle_t * p, p4est_quadrant_t * q);

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
static void
init_fn(p4est_t * p4est, p4est_topidx_t which_tree,
        p4est_quadrant_t * quadrant);
static int
refine_fn(p4est_t * p4est, p4est_topidx_t which_tree,
          p4est_quadrant_t * quadrant);
static int
                    coarsen_fn(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q[]);

/* prints all x,y data and pointers */
static void
                    viter_fn(p4est_iter_volume_info_t * info, void *Data);
/* frees all arrays from particles */
static void
                    destroy_fn(p4est_iter_volume_info_t * info, void *Data);
/* prints all particles in a xyz manner */
static void
                    print_fn(p4est_iter_volume_info_t * info, void *user_data);
/* moves all particles according to a given velocity field */
static void
                    update_x_fn(p4est_iter_volume_info_t * info, void *user_data);

/* returns true if particle lies in quad */
int
                    pchase_particle_lies_in_quad(const pchase_particle_t * p, p4est_quadrant_t * q);

/* returns true if particle lies in world */
int
pchase_particle_lies_in_world(pchase_world_t * W, const pchase_particle_t * p);

/* move particles from parent to children or vice versa */
static void
replace_fn(p4est_t * p4est, p4est_topidx_t which_tree,
           int num_outgoing, p4est_quadrant_t * outgoing[],
           int num_incoming, p4est_quadrant_t * incoming[]);

/* destroys all allocated data */
        void
                            pchase_world_destroy(pchase_world_t * W);
#endif
