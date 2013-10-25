#ifndef PCHASE_PARTICLE_H
#define PCHASE_PARTICLE_H

#include "p4est_extended.h"

/* this is the most basic particle flying through pchase */
typedef struct {
        double              x[DIM];
}
                    pchase_particle_t;

/** this holds all information a quadrant needs to know about our particles
 * nParticles	the number of all particles in this quadrant
 * p[25] 	place for 25 particles - more than 5 particles will not be
 * 		placed on purpose. This could happen only on particle movement
 * 
 *              Update: it never happened to exceed this bound, even not  in
 *              "big" simulations with 10^20, since the adaptive refinement
 *              (and accompanied by that the 2:1 relationship of the octants) 
 *              secures a partition fine enough to hold all particles.
 * 
 *              [of course there is an upper bound, where this won't work 
 *              any more - the finest possible refinement level 
 *              P4EST_QMAXLEVEL is actually set to 29
 *              => smalles cube has length 1, and biggest 2^29]
 */
typedef struct {
        int                 nParticles;
        pchase_particle_t   p[25];
}
                    pchase_quadrant_data_t;

#endif
