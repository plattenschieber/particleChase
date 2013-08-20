#define FILENAME_MAX 50

#include <p4est_bits.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pchase_world.h"
#include "pchase_particle.h"

typedef struct {
        MPI_Comm            mpicomm;
        int                 mpisize;
        int                 mpirank;
}
                    mpi_context_t;

int
main(int argc, char **argv)
{
        int                 mpiret;
        mpi_context_t       mpi_context, *mpi = &mpi_context;
        p4est_t            *p4est;
        p4est_connectivity_t *connectivity;
        pchase_world_t     *W;

        /* initialize MPI and p4est internals */
        mpiret = MPI_Init(&argc, &argv);
        SC_CHECK_MPI(mpiret);
        mpi->mpicomm = MPI_COMM_WORLD;  /* your favourite comm here */
        mpiret = MPI_Comm_size(mpi->mpicomm, &mpi->mpisize);
        SC_CHECK_MPI(mpiret);
        mpiret = MPI_Comm_rank(mpi->mpicomm, &mpi->mpirank);
        SC_CHECK_MPI(mpiret);

        /* Sets global program identifiers (e.g. the MPIrank) and some flags */
        sc_init(mpi->mpicomm, 1, 1, NULL, SC_LP_ALWAYS);
        /* Registers p4est with the SC Library and sets the logging behavior */
        p4est_init(NULL, SC_LP_DEFAULT);
        /* build up the world */
        W = pchase_world_init(p4est);

        /* store connectivity for a unitsquare */
        connectivity = p4est_connectivity_new_unitsquare();
        /* build uniform tree and get space for 25 particles each */
        p4est = p4est_new_ext(mpi->mpicomm, connectivity, 0, 1, 1,
                          sizeof(pchase_quadrant_data_t), W->init_fn, NULL);
        p4est_vtk_write_file(p4est, NULL, "pchase_new");

        /* don't forget to assign newly allocated p4est to the world */
        W->p4est = p4est;

        /* add one particle to the world */
        //pchase_world_insert_particle(W, pchase_world_random_particle(W));
        pchase_particle_t  *p = P4EST_ALLOC(pchase_particle_t, 1);
#ifdef DEBUG
        p->ID = 666;
#endif
        p->x[0] = 0.3;
        p->x[1] = 0.3;
        pchase_particle_t  *p2 = P4EST_ALLOC(pchase_particle_t, 1);
#ifdef DEBUG
        p2->ID = 666;
#endif
        p2->x[0] = 0.4;
        p2->x[1] = 0.4;
        pchase_particle_t  *p3 = P4EST_ALLOC(pchase_particle_t, 1);
#ifdef DEBUG
        p3->ID = 666;
#endif
        p3->x[0] = 0.91;
        p3->x[1] = 0.5;
        pchase_particle_t  *p4 = P4EST_ALLOC(pchase_particle_t, 1);
#ifdef DEBUG
        p4->ID = 666;
#endif
        p4->x[0] = 0.8;
        p4->x[1] = 0.5;
        pchase_particle_t  *p5 = P4EST_ALLOC(pchase_particle_t, 1);
#ifdef DEBUG
        p5->ID = 666;
#endif
        p5->x[0] = 0.7;
        p5->x[1] = 0.5;
        pchase_world_insert_particle(W, p);
        pchase_world_insert_particle(W, p2);
        pchase_world_insert_particle(W, p3);
        pchase_world_insert_particle(W, p4);
        pchase_world_insert_particle(W, p5);

        /* let this particle move */
        pchase_world_simulate(W);

#ifdef DEBUG
        /* print out all quadrants */
        p4est_iterate(W->p4est, NULL, NULL, W->viter_fn, NULL, NULL);
#endif

        /* destroy all particles, p4est and its connectivity structure */
        p4est_iterate(W->p4est, NULL, NULL, W->destroy_fn, NULL, NULL);
        p4est_destroy(p4est);
        p4est_connectivity_destroy(connectivity);

        /* clean up and exit */
        sc_finalize();
        mpiret = MPI_Finalize();
        SC_CHECK_MPI(mpiret);

        return 0;
}
