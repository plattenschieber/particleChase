#include <stdio.h>
#include <stdlib.h>
#include "particle.h"
#include "world.h"

int main() {
	PARTICLE p1;
   	PARTICLE p2;
	WORLD *W = world_init(NULL);

	/* insert first particle */
	p1.x[0] = 1.3; p1.x[1] = 2.3;
	p1.ID = 2;
	list_insert(W->cells, NULL, &p1);

	/* insert second particle on first place */
	p2.x[0] = 1.5; p2.x[1] = 4.8;
	p2.ID = 1;
	list_insert(W->cells,W->cells->first, &p2);

	/* print particles before and after the simulation */ 
	print_particles(W->cells);
	simulate(W);
	print_particles(W->cells);

	return 0;
}
