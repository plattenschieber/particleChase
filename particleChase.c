#define FILENAME_MAX 50

#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "list.h"
#include "world.h"

typedef struct
{
  MPI_Comm            mpicomm;
  int                 mpisize;
  int                 mpirank;
}
mpi_context_t;

	++PARTICLE_COUNT;
/* refines the cell only for first particle */
static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
	if (((PARTICLE *) (quadrant->p.user_data))->ID == 0)
	{
		return 1;
	}

	return 0;
}

int main(int argc, char **argv) {

	int                 mpiret;
	mpi_context_t       mpi_context, *mpi = &mpi_context;
	p4est_t            *p4est;
	p4est_connectivity_t *connectivity;
	p4est_refine_t      refine_fn;
	p4est_coarsen_t     coarsen_fn;

	/* initialize MPI and p4est internals */
	mpiret = MPI_Init (&argc, &argv);
	SC_CHECK_MPI (mpiret);
	mpi->mpicomm = MPI_COMM_WORLD;        /* your favourite comm here */
	mpiret = MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
	SC_CHECK_MPI (mpiret);
	mpiret = MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
	SC_CHECK_MPI (mpiret);
	sc_init (mpi->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
	p4est_init (NULL, SC_LP_DEFAULT);

	connectivity = p4est_connectivity_new_unitsquare ();
	

	/* argument handling */
	char file1[FILENAME_MAX];
	char file2[FILENAME_MAX];
	if (argc == 1) {
		strcpy(file1,"unit.particles");
		strcpy(file2,"unit.parameter");
#ifdef DEBUG
		printf("No parameter file, nor particles given - using 'unit.parameter' and 'unit.particles' instead\n");
#endif
	}
	else if (argc == 2) { 
#ifdef DEBUG
		printf("argv[1] = %s\n",argv[1]);
#endif
		strcpy(file1,argv[1]);
		strcpy(file2, "unit.parameter");
#ifdef DEBUG
		printf("No parameter file given, using 'unit.parameter' instead\n");
#endif
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
	print_particles(W->particles);
	world_simulate(W);
	print_particles(W->particles);

	return 0;
}
