#define FILENAME_MAX 50

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "list.h"
#include "world.h"

int main(int argc, char **argv) {
	/* argument handling */
	char file1[FILENAME_MAX];
	char file2[FILENAME_MAX];
	if (argc == 1) {
		strcpy(file1,"unit.particles");
		strcpy(file2,"unit.parameter");
		printf("No parameter file, nor particles given - using 'unit.parameter' and 'unit.particles' instead\n");
	}
	else if (argc == 2) { 
#ifdef DEBUG
#endif
		strcpy(file1,argv[1]);
		strcpy(file2, "unit.parameter");
		printf("No particles file given, using 'unit.particles' instead\n");
	}
	else if (argc == 3) { 
		strcpy(file1,argv[1]);
		strcpy(file2,argv[2]);
	}
	else 
		printf("Too many arguments.\n");

	/* open file handling */
	FILE *parameter, *particles; 
	if ((particles = fopen(file1,"r"))==NULL) {
		fprintf(stderr, "Could not open %s. Please Check your file.\n", file2);
		exit(EXIT_FAILURE); 
	}
	if ((parameter = fopen(file2,"r"))==NULL) {
		fprintf(stderr, "Could not open %s. Please Check your file.\n", file1);
		exit(EXIT_FAILURE); 
	}

	/* initialize the world and close files */
	WORLD *W = world_init(parameter, particles);
	fclose(parameter);
	fclose(particles);
	/* print particles before and after the simulation */ 
	print_particles(W->cells);
	world_simulate(W);
	print_particles(W->cells);

	return 0;
}
