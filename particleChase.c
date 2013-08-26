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
        W = pchase_world_init();

        /* store connectivity for a unitsquare */
        connectivity = p4est_connectivity_new_unitsquare();
        /* build uniform tree and get space for 25 particles each */
        p4est = p4est_new_ext(mpi->mpicomm, connectivity, 0, 2, 1,
                          sizeof(pchase_quadrant_data_t), W->init_fn, NULL);

        /* initialize everything depending on p4est */
        pchase_world_init_p4est(W, p4est);

        if (W->p4est->mpirank == 0) {
                pchase_particle_t  *p = P4EST_ALLOC(pchase_particle_t, 1);
                p->ID = 1;
                p->x[0] = 0.3;
                p->x[1] = 0.3;
                sc_list_append(W->particle_push_list, p);
                W->n_particles++;
                pchase_particle_t  *p2 = P4EST_ALLOC(pchase_particle_t, 1);
                p2->ID = 2;
                p2->x[0] = 0.4;
                p2->x[1] = 0.4;
                sc_list_append(W->particle_push_list, p2);
                W->n_particles++;
                pchase_particle_t  *p3 = P4EST_ALLOC(pchase_particle_t, 1);
                p3->ID = 3;
                p3->x[0] = 0.9;
                p3->x[1] = 0.5;
                sc_list_append(W->particle_push_list, p3);
                W->n_particles++;
                pchase_particle_t  *p4 = P4EST_ALLOC(pchase_particle_t, 1);
                p4->ID = 4;
                p4->x[0] = 0.4;
                p4->x[1] = 0.9;
                sc_list_append(W->particle_push_list, p4);
                W->n_particles++;
                pchase_particle_t  *p5 = P4EST_ALLOC(pchase_particle_t, 1);
                p5->ID = 5;
                p5->x[0] = 0.8;
                p5->x[1] = 0.6;
                sc_list_append(W->particle_push_list, p5);
                W->n_particles++;
        }
        /* this has to be done for each proc */
        pchase_world_insert_particles(W);

        /* /1* add a bunch of particles to the world *1/ */
        /* for (i = 0; i < 400; i++) { */
        /* pchase_particle_t  *p = P4EST_ALLOC(pchase_particle_t, 1); */
        /* #ifdef DEBUG */
        /* p->ID = i; */
        /* #endif */
        /* p->x[0] = i * 0.0001 + 0.55; */
        /* p->x[1] = 0.5; */
        /* sc_list_append(W->particle_push_list, p); */
        /* W->n_particles++; */
        /* } */
        /* pchase_world_insert_particles(W); */

        /* let this particle move */
        pchase_world_simulate(W);

#ifdef DEBUG
        /* print out all quadrants */
        p4est_iterate(W->p4est, NULL, NULL, W->viter_fn, NULL, NULL);
#endif

        /* destroy all particles, p4est and its connectivity structure */
        pchase_world_destroy(W);
        p4est_destroy(p4est);
        p4est_connectivity_destroy(connectivity);

        /* clean up and exit */
        sc_finalize();
        mpiret = MPI_Finalize();
        SC_CHECK_MPI(mpiret);

        return 0;
}
