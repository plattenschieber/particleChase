#ifndef PCHASE_PARTICLE_H
#define PCHASE_PARTICLE_H
#include "world.h"

typedef struct {
#ifdef DEBUG /* needed only to follow the particle */
	int ID;
#endif
	double x[DIM];
} pchase_particle_t;

typedef struct {
  int nParticles;
  pchase_particle_t p[25];
} pchase_quadrant_data_t;

#endif

