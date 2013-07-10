#include <stdio.h>
#include <stdlib.h>
#include "list.h"
#include "world.h"

int main() {


	/* initialize the world and close files */
	WORLD *W = world_init(parameter, particles);
	/* print particles before and after the simulation */ 
	print_particles(W->cells);
	simulate(W);
	print_particles(W->cells);

	return 0;
}
