#ifndef PCHASE_PARTICLE_H
#define PCHASE_PARTICLE_H

#include "p4est_extended.h"

/* this is the most basic particle flying through pchase */
typedef struct {
#ifdef DEBUG                    /* needed only to follow the particle */
        int                 ID;
#endif
        double              x[DIM];
}
                    pchase_particle_t;

/** this holds all information a quadrant needs to know about our particles
 * nParticles	the number of all particles in this quadrant
 * p[25] 	place for 25 particles - more than 5 particles will not be
 * 		placed on purpose. This could happen only on particle movement
 */
typedef struct {
        int                 nParticles;
        pchase_particle_t  p[25];
}
                    pchase_quadrant_data_t;

#endif
