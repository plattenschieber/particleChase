#define FILENAME_MAX 50

#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pchase_particle.h"
#include "pchase_world.h"

typedef struct
{
	MPI_Comm            mpicomm;
	int                 mpisize;
	int                 mpirank;
}
mpi_context_t;

/* every quadrant will get one particle with ID -1 */
static void
init_fn (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
	quadrant->p.user_data = NULL;
}

static int
refine_fn (p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
	if (0)
	  return 1;
	else 
	  return 0;
}

int main(int argc, char **argv) {

  int                 mpiret;
  mpi_context_t       mpi_context, *mpi = &mpi_context;
  p4est_t            *p4est;
  p4est_connectivity_t *connectivity;
  p4est_refine_t      refine_fn;

  /* initialize MPI and p4est internals */
  mpiret = MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpi->mpicomm = MPI_COMM_WORLD;        /* your favourite comm here */
  mpiret = MPI_Comm_size (mpi->mpicomm, &mpi->mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = MPI_Comm_rank (mpi->mpicomm, &mpi->mpirank);
  SC_CHECK_MPI (mpiret);

  /* Sets the global program identifier (e.g. the MPI rank) and some flags */
  sc_init (mpi->mpicomm, 1, 1, NULL, SC_LP_ALWAYS); 
  /* Registers p4est with the SC Library and sets the logging behavior */
  p4est_init (NULL, SC_LP_DEFAULT);

  /* store connectivity for a unitsquare */
  connectivity = p4est_connectivity_new_unitsquare ();

  /* build up tree and require quadrant space for 25 particles each  */
  p4est = p4est_new_ext (mpi->mpicomm, connectivity, 0, 0, 0,
	  sizeof (pchase_quadrant_data_t), init_fn, NULL);

  p4est_vtk_write_file (p4est, NULL, "pchase_new");

  /* destroy the p4est and its connectivity structure */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (connectivity);

  /* clean up and exit */
  sc_finalize ();

  mpiret = MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
