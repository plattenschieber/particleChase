#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "list.h"
#include "world.h"

int main(int argc, char **argv) {
	/* argument handling */
	char file1[100];
	char file2[100];
	if (argc == 1) {
		strcpy(file1,"unit.parameter");
		strcpy(file2,"unit.particles");
		printf("No parameter file, nor particles given - using 'unit.parameter' and 'unit.particles' instead\n");
	}
	else if (argc == 2) { 
		strcpy(file1,argv[2]);
		strcpy(file2, "unit.parameter");
		printf("No parameter file given, using 'unit.parameter' instead\n");
	}
	else if (argc == 3) { 
		strcpy(file1,argv[2]);
		strcpy(file2,argv[3]);
	}
	else 
		printf("Too many arguments.\n");


	/* initialize the world and close files */
	WORLD *W = world_init(parameter, particles);
	/* print particles before and after the simulation */ 
	print_particles(W->cells);
	simulate(W);
	print_particles(W->cells);

	return 0;
}
