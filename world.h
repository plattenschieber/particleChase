#ifndef WORLD_H
#define WORLD_H
#include "particle.h"

/* a list of cells with particle lists */
LIST *cells;
// current time */
double world_t;
// timestep */
double world_delta_t;
// end time */
double world_t_end;
// number of overall particles */
unsigned int n_particles;
/* start the simulation */
void simulate(LIST *L);
/* update the positions */
void update_x();

#endif