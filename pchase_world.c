#include "pchase_world.h"

static FILE        *pchase_output;

pchase_world_t     *
pchase_world_init(p4est_t * p4est)
{
        int                 i;
        /* get some place for the pchase_world */
        pchase_world_t     *W;
        if (!(W = malloc(sizeof(pchase_world_t)))) {
                fprintf(stderr, "Could not allocate memory for the World\n");
                exit(EXIT_FAILURE);
        }
        /* set all parameters */
        W->t = 0.0;
        W->delta_t = 0.001;
        W->t_end = 0.500;
        W->n_particles = 0;
        W->step = 0;

        /* set callback functions for various applications */
        W->init_fn = init_fn;
        W->coarsen_fn = coarsen_fn;
        W->refine_fn = refine_fn;
        W->search_fn = search_fn;
        W->replace_fn = replace_fn;
        W->viter_fn = viter_fn;
        W->print_fn = print_fn;
        W->update_x_fn = update_x_fn;

        /* take care of all particles which left the cell after update_x */
        W->particle_push_list = sc_list_new(NULL);
        /* set world length */
        for (i = 0; i < DIM; i++)
                W->length[i] = 1.0;
        /* reset seed */
        srand(time(NULL));
        return W;
}

void
pchase_world_init_p4est(pchase_world_t * W, p4est_t * p4est)
{
        int                 i;
        /* don't forget to assign newly allocated p4est to the world */
        W->p4est = p4est;
        /* initialize particle send list */
        W->particles_to = sc_array_new_size(sizeof(sc_list_t *), W->p4est->mpisize);
#ifdef PRINTXYZ
        pchase_output = fopen("pchase_particles.xyz", "w");
#elif defined(PRINTGNUPLOT)
        /* generate extension for proc i and open its file */
        char                plotFileName[50] = "pchase_particles_";
        char                fileNumber[10];
        sprintf(fileNumber, "%03d", W->p4est->mpirank);
        strcat(plotFileName, fileNumber);
        strcat(plotFileName, ".plot");
        pchase_output = fopen(plotFileName, "w");
#endif
        /* reserve some space for send lists */
        for (i = 0; i < W->p4est->mpisize; i++)
                *((sc_list_t **) sc_array_index_int(W->particles_to, i)) = sc_list_new(NULL);
}

void
pchase_world_simulate(pchase_world_t * W)
{
        int                 i = 0;
        FILE               *vtk_timeseries = fopen("pchase_particle_simulation.pvd", "w");

        fprintf(vtk_timeseries, "<?xml version=\"1.0\"?>\n");
        fprintf(vtk_timeseries, "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
        fprintf(vtk_timeseries, "      <Collection>\n");

        char                fileName1[100] = "pchase_with_particles_";
        char                fileName[100] = "";
        char                VTKData1[100] = "            <DataSet timestep=\"";
        char                VTKData2[100] = "\" file=\"";
        char                VTKData3[100] = ".pvtu\"/>\n";
        char                VTKData[200] = "";
        char                fileNumber[10];

        printf("[pchase %i simulate] start new simulation\n", W->p4est->mpirank);
        /* p4est_iterate(W->p4est, NULL, W, W->print_fn, NULL, NULL); */
        /* simulate until the end has come */
        while (W->t <= W->t_end) {
#ifdef PRINTXYZ
                if (W->step % 1000 == 0) {
                        fprintf(pchase_output, "%i\n", W->n_particles);
                        fprintf(pchase_output, "Time: %f\n", W->t);
                }
#endif
#ifdef DEBUG
                printf("[pchase %i simulate] into update_x\n", W->p4est->mpirank);
#endif
                /* update the position of all particles on all quads */
                p4est_iterate(W->p4est, NULL, W, W->update_x_fn, NULL, NULL);
#ifdef DEBUG
                printf("[pchase %i simulate] update_x done - starting insertion of particle\n", W->p4est->mpirank);
#endif
                /* insert all particles which left into their new quad */
                pchase_world_insert_particles(W);
                /* insert all communicated particles into the world */
                pchase_world_insert_particles(W);
                /* refine every quad containing more than 5 particles */
                p4est_refine_ext(W->p4est, 0, -1, W->refine_fn, W->init_fn, W->replace_fn);
                /* coarsen non recursively */
                p4est_coarsen_ext(W->p4est, 0, W->coarsen_fn, W->init_fn, W->replace_fn);
                /* balancing the tree */
                p4est_balance_ext(W->p4est, P4EST_CONNECT_FULL, W->init_fn, W->replace_fn);
                /* partition quadrants evenly to all procs */ 
                p4est_partition_ext(W->p4est, 1, NULL);

                /*
                 * convert current step to filename and write VTK
                 * Data entry
                 */
                if (W->step % 50 == 0) {
                        sprintf(fileNumber, "%03d", i);
                        strcat(VTKData, VTKData1);
                        strcat(VTKData, fileNumber);
                        strcat(VTKData, VTKData2);
                        strcat(fileName, fileName1);
                        strcat(fileName, fileNumber);
                        strcat(VTKData, fileName);
                        strcat(VTKData, VTKData3);
                        fprintf(vtk_timeseries, "%s", VTKData);
                        p4est_vtk_write_file(W->p4est, NULL, fileName);

                        i++;
                        /* reset filenames */
                        VTKData[0] = '\0';
                        fileName[0] = '\0';
                }
                W->t += W->delta_t;
                W->step++;
        }
        printf("[pchase %i simulate] simulation over\n", W->p4est->mpirank);
        /* p4est_iterate(W->p4est, NULL, W, W->print_fn, NULL, NULL); */
        fprintf(vtk_timeseries, "      </Collection>\n");
        fprintf(vtk_timeseries, "</VTKFile>\n");
        fclose(vtk_timeseries);
        fclose(pchase_output);
}

pchase_particle_t  *
pchase_world_random_particle(pchase_world_t * W)
{
        int                 i;
        double x,y;
        pchase_particle_t  *p = P4EST_ALLOC(pchase_particle_t, 1);

        do {
                for (i = 0; i < DIM; i++)
                        p->x[i] = W->length[i] * rand() / (RAND_MAX + 1.);
                x = p->x[0]-0.5;
                y = p->x[1]-0.5;
        } while (x*x+y*y >= 0.25 || x*x+y*y<0.05);
        
#ifdef DEBUG
        printf("[pchase %i randPart] ", W->p4est->mpirank);
        for (i = 0; i < DIM; i++)
                printf("p->x[%d]=%lf,\t", i, p->x[i]);
        printf("\n");
#endif
        return p;
}

void
pchase_world_insert_particles(pchase_world_t * W)
{
        p4est_quadrant_t   *miniQuad;
        pchase_particle_t  *p;
        sc_array_t         *points;
        int                 owner, i;
        int                 mpiret;
        int                *receivers, num_receivers;
        int                *senders, num_senders;

        /* get place for as many points as particles waiting to be inserted */
        points = sc_array_new_size(sizeof(p4est_quadrant_t), W->particle_push_list->elem_count);
        /* set the particle iterator */
        sc_link_t          *particle_it = W->particle_push_list->first;

        /* calculate a miniQuad for every particle in particle_push_list */
        for (i = 0; i < points->elem_count; i++, particle_it = particle_it->next) {
                /* resolve to be filled miniQuad and particle */
                miniQuad = (p4est_quadrant_t *) sc_array_index(points, i);
                p = (pchase_particle_t *) particle_it->data;

                /* create miniQuad that is enclosing the given particle p */
                pchase_translate_particle_to_p4est(W, p, miniQuad);

#ifdef DEBUG
                printf("[pchase %i insertPart] Translated P(%lf,%lf) to miniQuad (0x%08X,0x%08X)\n",
                       W->p4est->mpirank, p->x[0], p->x[1], miniQuad->x, miniQuad->y);
#endif
        }

        /*
         * find all quads enclosing the miniQuads and save their data in the
         * miniQuads' piggy3
         */
        p4est_search(W->p4est, W->search_fn, points);
#ifdef DEBUG
        printf("[pchase %i insertPart] %lld particles in push list\n", W->p4est->mpirank, (long long)W->particle_push_list->elem_count);
#endif
        /* move all particles either into their enclQuads or to another proc */
        for (i = 0; i < points->elem_count; i++) {
                /* resolve miniQuad and associated particle */
                miniQuad = (p4est_quadrant_t *) sc_array_index(points, i);

                p = sc_list_pop(W->particle_push_list);

                /* if the miniQuad is not flagged, it's lying on this proc */
                if (miniQuad->p.piggy3.local_num != -1) {
                        /* extract enclosing quad from miniQuad.piggy3 */
                        p4est_tree_t       *enclQuadTree = p4est_tree_array_index(W->p4est->trees, miniQuad->p.piggy3.which_tree);
                        p4est_quadrant_t   *enclQuad = p4est_quadrant_array_index(&enclQuadTree->quadrants, miniQuad->p.piggy3.local_num);
                        pchase_quadrant_data_t *enclQuadData = enclQuad->p.user_data;

                        if (enclQuadData->nParticles < 25) {
                                /*
                                 * copy particle data into quad and update
                                 * particle counter
                                 */
                                enclQuadData->p[enclQuadData->nParticles] = *p;
                                enclQuadData->nParticles++;
#ifdef DEBUG
                                printf("[pchase %i insertPart] P(%lf,%lf) inserted into enclQuad(0x%08X,0x%08X) (had already %i particles - now %i)\n",
                                       W->p4est->mpirank, p->x[0], p->x[1], enclQuad->x, enclQuad->y, enclQuadData->nParticles - 1, enclQuadData->nParticles);
#endif
                                P4EST_FREE(p);
                        }
                        /* too many particles in quad */
                        else {
                                printf("[pchase %i insertPart] Too many (%i) particles ",
                                W->p4est->mpirank, enclQuadData->nParticles);
                                if (0) {
                                        printf("- refining enclQuad and inserting particle afterwards\n");
                                        p4est_refine_ext(W->p4est, 0, -1, W->refine_fn, W->init_fn, W->replace_fn);
                                        /*
                                         * this can be done even in this
                                         * special case, since we are
                                         * inserting as long as this list is
                                         * not empty
                                         */
                                        sc_list_append(W->particle_push_list, p);
                                } else {
                                        printf("- we have to dissmiss this particle\n");
                                        P4EST_FREE(p);
                                }
                        }
                }
                /* particle lies on another proc */
                else {
                        /* resolving particles' owner */
                        owner = p4est_comm_find_owner(W->p4est, miniQuad->p.piggy3.which_tree, miniQuad, W->p4est->mpirank);
                        /* moving particle into sent list for proc 'owner' */
                        sc_list_t          *tmp = *((sc_list_t **) sc_array_index(W->particles_to, owner));
#ifdef DEBUG
                        printf("[pchase %i insertPart] particle(%lf,%lf) pushed to send list for proc %i (had already %lld particles - ",
                               W->p4est->mpirank, p->x[0], p->x[1], owner, (long long)tmp->elem_count);
#endif
                        sc_list_append(tmp, p);
#ifdef DEBUG
                        printf("now %lld)\n", (long long)tmp->elem_count);
#endif
                }
        }

        /* WE ARE DONE WITH ALL LOCAL STUFF */

        /* reset push list - we need this for saving incomming particles */
        sc_list_reset(W->particle_push_list);

        /*
         * get enough space for receivers and senders array - this may be not
         * memory optimal, but it's the fastest solution
         */
        receivers = P4EST_ALLOC(int, W->p4est->mpisize);
        senders = P4EST_ALLOC(int, W->p4est->mpisize);
        /* resolve receiver count */
        for (i = 0, num_receivers = 0; i < W->p4est->mpisize; i++) {
                sc_list_t          *tmp = *((sc_list_t **) sc_array_index(W->particles_to, i));
                /* add i to receiver list and update receiver count */
                if (tmp->elem_count)
                        receivers[num_receivers++] = i;
        }
#ifdef DEBUG
        printf("[pchase %i sending] total num_receivers %i \n", W->p4est->mpirank, num_receivers);
        /* printing all particles to be sent */
        for (i = 0; i < W->particles_to->elem_count; i++) {
                sc_list_t          *tmpList = *((sc_list_t **) sc_array_index(W->particles_to, i));
                if (tmpList->elem_count > 0) {
                        printf("[pchase %i sending] %lld particles to proc %i: ", W->p4est->mpirank, (long long)tmpList->elem_count, i);
                        sc_link_t          *tmpLink = tmpList->first;
                        pchase_particle_t  *tmpParticle;
                        while (tmpLink != NULL) {
                                tmpParticle = tmpLink->data;
                                printf("P(%lf,%lf) ", tmpParticle->x[0], tmpParticle->x[1]);
                                tmpLink = tmpLink->next;
                        }
                        printf("\n");
                }
        }
#endif

        /* only the first num_receivers receivers are notified */
        mpiret = sc_notify(receivers, num_receivers,
                           senders, &num_senders, W->p4est->mpicomm);
        SC_CHECK_MPI(mpiret);

#ifdef DEBUG
        printf("[pchase %i insertPart] sc_notify done num_receivers %i, num_senders %i\n", W->p4est->mpirank, num_receivers, num_senders);
#endif

        /* handle all mpi send/recv status data */
        MPI_Request        *send_request = P4EST_ALLOC(MPI_Request, W->p4est->mpisize);
        MPI_Status         *recv_status = P4EST_ALLOC(MPI_Status, W->p4est->mpisize);
        /* setup send/recv buffers */
        double            **recv_buf = P4EST_ALLOC(double *, num_senders);
        double            **send_buf = P4EST_ALLOC(double *, num_receivers);
        int                 recv_count = 0, recv_length, flag, j;

        /* send all particles to their belonging procs */
        for (i = 0; i < num_receivers; i++) {
                /* resolve particle list for proc i */
                sc_list_t          *tmpList = *((sc_list_t **) sc_array_index(W->particles_to, receivers[i]));
                pchase_particle_t  *tmpParticle;
                int                 send_count = 0;

                /* get space for the particles to be sent */
                send_buf[i] = P4EST_ALLOC(double, 2 * tmpList->elem_count);

                /*
                 * copy all particles into the send buffer and remove them
                 * from this proc
                 */
                while (tmpList->first != NULL) {
                        tmpParticle = sc_list_pop(tmpList);
                        send_buf[i][2 * send_count] = tmpParticle->x[0];
                        send_buf[i][2 * send_count + 1] = tmpParticle->x[1];
                        /* free particle */
                        P4EST_FREE(tmpParticle);
                        /* update particle counter */
                        send_count++;
                }

#ifdef DEBUG
                /* print send buffer */
                for (j = 0; j < send_count; j++) {
                        printf("[pchase %i sending] P(%lf,%lf)\n",
                               W->p4est->mpirank, send_buf[i][2 * j], send_buf[i][2 * j + 1]);
                }
#endif
                /* send particles to right owner */
                mpiret = MPI_Isend(send_buf[i], 2 * send_count, MPI_DOUBLE,
                                   receivers[i], 13,
                                   W->p4est->mpicomm, &send_request[i]);
                SC_CHECK_MPI(mpiret);
        }

        recv_count = 0;
        /* check for messages until all arrived */
        while (recv_count < num_senders) {
                /* probe if any of the sender has already sent his message */
                for (i = 0; i < num_senders; i++) {
                        MPI_Iprobe(senders[i], MPI_ANY_TAG, W->p4est->mpicomm,
                                   &flag, &recv_status[i]);
                        if (flag) {
                                /* resolve number of particles receiving */
                                MPI_Get_count(&recv_status[i], MPI_DOUBLE, &recv_length);
#ifdef DEBUG
                                printf("[pchase %i receiving] %i particles from proc %i: ",
                                       W->p4est->mpirank, recv_length / 2, recv_status[i].MPI_SOURCE);
#endif
                                /* get space for the particles to be received */
                                recv_buf[recv_count] = P4EST_ALLOC(double, recv_length);
                                /*
                                 * receive a list with recv_length/2
                                 * particles
                                 */
                                mpiret = MPI_Recv(recv_buf[recv_count], recv_length, MPI_DOUBLE, recv_status[i].MPI_SOURCE,
                                                  recv_status[i].MPI_TAG, W->p4est->mpicomm, &recv_status[i]);
                                SC_CHECK_MPI(mpiret);

                                /*
                                 * insert all received particles into the
                                 * push list
                                 */
                                pchase_particle_t  *tmpParticle;
                                for (j = 0; j < recv_length / 2; j++) {
                                        /*
                                         * retrieve all particle details from
                                         * recv_buf
                                         */
                                        tmpParticle = P4EST_ALLOC(pchase_particle_t, 1);
                                        tmpParticle->x[0] = recv_buf[recv_count][2 * j];
                                        tmpParticle->x[1] = recv_buf[recv_count][2 * j + 1];

#ifdef DEBUG
                                        printf("Particle(%lf,%lf) ",
                                               tmpParticle->x[0], tmpParticle->x[1]);
#endif
                                        /*
                                         * push received particle to push
                                         * list and update world counter
                                         */
                                        sc_list_append(W->particle_push_list, tmpParticle);
                                }
#ifdef DEBUG
                                printf("\n");
#endif
                                /* we received another particle list */
                                recv_count++;
                        }
                }
        }

        /*
         * wait for all procs to finish sending (recieve is blocking, so no
         * need to wait here)
         */
        if (num_receivers > 0) {
                /* wait for receivers to collect all messages */
                mpiret = MPI_Waitall(num_receivers, send_request, MPI_STATUSES_IGNORE);
                SC_CHECK_MPI(mpiret);
        }
        /* free proc lists */
        P4EST_FREE(receivers);
        P4EST_FREE(senders);
        /* free mpi handles */
        P4EST_FREE(send_request);
        P4EST_FREE(recv_status);
        /* get rid of receive and sendbuffer */ 
        for(i=0; i<num_senders; i++)
                P4EST_FREE(recv_buf[i]);
        P4EST_FREE(recv_buf);
        for(i=0; i<num_receivers; i++)
                P4EST_FREE(send_buf[i]);
        P4EST_FREE(send_buf);
        /* free all miniQuads */ 
        sc_array_destroy(points);
}

void
pchase_translate_particle_to_p4est(pchase_world_t * W, const pchase_particle_t * p, p4est_quadrant_t * q)
{
        double              quadrant_length = (double)(1 << P4EST_QMAXLEVEL);

        /*
         * normalize particle and transform it to p4est world length by
         * truncating its position on the quadrant grid and MULT by 2, to get
         * its real position in p4est
         */
        q->x = (p4est_qcoord_t) (p->x[0] / W->length[0] * quadrant_length) << 1;
        q->y = (p4est_qcoord_t) (p->x[1] / W->length[1] * quadrant_length) << 1;
#if DIM == 3
        q->z = (p4est_qcoord_t) (p->x[2] / W->length[2] * quadrant_length) << 1;
#endif
        q->level = P4EST_MAXLEVEL;
}

void
pchase_world_print_particlesXYZ(pchase_world_t * W)
{
        printf("PRINT XYZ - NOT IMPLEMENTED YET");
}

int
pchase_quadrant_is_in_proc(p4est_quadrant_t * q)
{
        return -1;
}

static int
search_fn(p4est_t * p4est, p4est_topidx_t which_tree,
          p4est_quadrant_t * quadrant, int is_leaf, void *point)
{
        p4est_quadrant_t   *miniQuad;
        int                 quadrant_length;

        quadrant_length = P4EST_QUADRANT_LEN(quadrant->level);
        miniQuad = (p4est_quadrant_t *) point;

        /* mini quad must lie entirely in quadrant */
        if (miniQuad->x >= quadrant->x && miniQuad->x < quadrant->x + quadrant_length &&
            miniQuad->y >= quadrant->y && miniQuad->y < quadrant->y + quadrant_length) {
                /* remember current tree (needed to identify proc) */
                miniQuad->p.piggy3.which_tree = which_tree;
                if (is_leaf) {
                        /*
                         * miniQuad lies on this proc, local num holds info
                         */
                        miniQuad->p.piggy3.local_num = quadrant->p.piggy3.local_num;
#ifdef DEBUG
                        printf("[pchase %i search] found enclQuad[%d] on level: %d\n",
                               p4est->mpirank, miniQuad->p.piggy3.local_num, quadrant->level);
#endif
                } else
                        /* flag the quad - it's lying on another proc */
                        /* (or still not the leaf, but keep flagging) */
                        miniQuad->p.piggy3.local_num = -1;
                return 1;
        } else
                return 0;
}

static void
init_fn(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
#ifdef DEBUG
        printf("[pchase %i init_fn] quad(0x%08X,0x%08X) in tree: %lld\n", p4est->mpirank, quadrant->x, quadrant->y, (long long)which_tree);
#endif
        ((pchase_quadrant_data_t *) quadrant->p.user_data)->nParticles = 0;
}

static int
refine_fn(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
        if (((pchase_quadrant_data_t *) quadrant->p.user_data)->nParticles > 5)
                return 1;
        else
                return 0;
}

static int
coarsen_fn(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * q[])
{
        int                 i, sum;
        for (i = 0, sum = 0; i < P4EST_CHILDREN; i++)
                sum += ((pchase_quadrant_data_t *) q[i]->p.user_data)->nParticles;

        if (sum > 5)
                return 0;
        else
                return 1;
}

static void
viter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
        pchase_quadrant_data_t *quadData = (pchase_quadrant_data_t *) info->quad->p.user_data;
        if (quadData->nParticles > 0)
                printf("[pchase %i main iterate] quad(0x%08X,0x%08X) has %i particles \n",
                       info->p4est->mpirank, info->quad->x, info->quad->y, quadData->nParticles);
        return;
}
static void
print_fn(p4est_iter_volume_info_t * info, void *user_data)
{
        int                 i;
        pchase_quadrant_data_t *quadData = (pchase_quadrant_data_t *) info->quad->p.user_data;
        if (quadData->nParticles > 0) {
                printf("[pchase %i print] quad[%i](0x%08X,0x%08X)", info->p4est->mpirank, quadData->nParticles, info->quad->x, info->quad->y);
                for (i = 0; i < quadData->nParticles; i++)
                        printf("P(%lf,%lf) ", quadData->p[i].x[0], quadData->p[i].x[1]);
                printf("\n");
                fflush(stdout);
        }
}

void
pchase_world_velocity(pchase_world_t * W, pchase_particle_t * p)
{
        double              x, y, norm;
        x = -p->x[1] + 0.5;
        y = p->x[0] - 0.5;

        norm = sqrt(x * x + y * y);
        p->x[0] += W->delta_t * x / norm;
        p->x[1] += W->delta_t * y / norm;

        if (!pchase_particle_lies_in_world(W, p))
                sc_abort_collective("Particle is lying inside the World");
}

static void
update_x_fn(p4est_iter_volume_info_t * info, void *user_data)
{
        int                 i;
        pchase_quadrant_data_t *quadData = (pchase_quadrant_data_t *) info->quad->p.user_data;
        pchase_world_t     *W = (pchase_world_t *) user_data;

#ifdef DEBUG
        if (quadData->nParticles > 0)
                printf("[pchase %i updateX] quad(0x%08X,0x%08X) with %i P: ",
                       W->p4est->mpirank, info->quad->x, info->quad->y, quadData->nParticles);
#endif
        for (i = 0; i < quadData->nParticles; i++) {
#ifdef DEBUG
                /* print old position */
                printf("P(%lf,%lf)->",
                       quadData->p[i].x[0], quadData->p[i].x[1]);
#endif

                /* update particles' velocity */
                pchase_world_velocity(W, &quadData->p[i]);

#ifdef DEBUG
                /* and new position */
                printf("(%lf,%lf) ",
                       quadData->p[i].x[0], quadData->p[i].x[1]);
#endif
#ifdef PRINTXYZ
                fprintf(pchase_output, "H\t%lf\t%lf\t0\n", quadData->p[i]->x[0], quadData->p[i]->x[1]);
#elif defined(PRINTGNUPLOT)
                fprintf(pchase_output, "%lf\t%lf\n", quadData->p[i].x[0], quadData->p[i].x[1]);
#endif
                /* move particle if it has left the quad */
                if (!pchase_particle_lies_in_quad(&quadData->p[i], info->quad)) {
#ifdef DEBUG
                        printf(" (left quad!)");
#endif
                        /*
                         * pop particle and leave a pointer to it in
                         * particle_push_list
                         */
                        pchase_particle_t  *tmpParticle = P4EST_ALLOC(pchase_particle_t, 1);
                        *tmpParticle = quadData->p[i];
                        sc_list_append(W->particle_push_list, tmpParticle);

                        /*
                         * move last particle to i'th place to prevent a hole
                         */
                        quadData->p[i] = quadData->p[quadData->nParticles - 1];
                        /* set iterator accordingly */
                        i--;

                        /* update particle counter in quad array */
                        quadData->nParticles--;
                }
        }
#ifdef DEBUG
        if (quadData->nParticles > 0)
                printf("\n");
#endif
}

int
pchase_particle_lies_in_world(pchase_world_t * W, const pchase_particle_t * p)
{
        if (p->x[0] < 0 || p->x[0] >= W->length[0] ||
            p->x[1] < 0 || p->x[1] >= W->length[1]) {
                printf("P(%lf,%lf) left World\n", p->x[0], p->x[1]);
                return 0;
        } else
                return 1;
}

int
pchase_particle_lies_in_quad(const pchase_particle_t * p, p4est_quadrant_t * q)
{
        p4est_qcoord_t      quadrant_length, root_len;
        quadrant_length = P4EST_QUADRANT_LEN(q->level);
        root_len = P4EST_ROOT_LEN;

        if (p->x[0] * root_len < q->x || p->x[0] * root_len >= q->x + quadrant_length ||
            p->x[1] * root_len < q->y || p->x[1] * root_len >= q->y + quadrant_length) {
                return 0;
        } else
                return 1;
}
static void
replace_fn(p4est_t * p4est, p4est_topidx_t which_tree,
           int num_outgoing, p4est_quadrant_t * outgoing[],
           int num_incoming, p4est_quadrant_t * incoming[])
{
        p4est_quadrant_t  **fam;
        p4est_quadrant_t   *p;
        pchase_quadrant_data_t *quadData, *famJData;
        int                 i, j;

        /* refining quad -> spread data to its children */
        if (num_outgoing == 1) {
#ifdef DEBUG
                printf("[pchase %i replace] replacing quad by %i children\n", p4est->mpirank, P4EST_CHILDREN);
#endif
                /* set readable names */
                p = outgoing[0];
                quadData = (pchase_quadrant_data_t *) p->p.user_data;
                fam = incoming;
                for (i = 0; i < quadData->nParticles; i++)
                        for (j = 0; j < P4EST_CHILDREN; j++) {
                                if (pchase_particle_lies_in_quad(&quadData->p[i], fam[j])) {
                                        /*
                                         * move particle to this quad(fam[j])
                                         */
                                        famJData = (pchase_quadrant_data_t *) fam[j]->p.user_data;
                                        famJData->p[famJData->nParticles] = quadData->p[i];
                                        /*
                                         * move the last particle to i'th
                                         * place to prevent a hole
                                         */
                                        quadData->p[i] = quadData->p[quadData->nParticles - 1];
                                        /*
                                         * reset iterator and particle
                                         * counter
                                         */
                                        i--;
                                        famJData->nParticles++;
                                        quadData->nParticles--;
                                        break;
                                }
                        }
        }
        /*
         * coarsening children -> gather data from them and sent to parent
         */
        else {
#ifdef DEBUG
                printf("[pchase %i replace] replacing %i children by their parent: ", p4est->mpirank, P4EST_CHILDREN);
#endif
                /* set readable names */
                p = incoming[0];
                quadData = (pchase_quadrant_data_t *) p->p.user_data;
                fam = outgoing;
                for (j = 0; j < P4EST_CHILDREN; j++) {
                        famJData = (pchase_quadrant_data_t *) fam[j]->p.user_data;
#ifdef DEBUG
                        printf("quad[%i](0x%08X,0x%08X)[%i] ", j, fam[j]->x, fam[j]->y, famJData->nParticles);
#endif
                        for (i = 0; i < famJData->nParticles; i++) {
                                /* move particle to parent */
                                quadData->p[quadData->nParticles] = famJData->p[i];
                                /*
                                 * move last particle in array to originated
                                 * hole
                                 */
                                famJData->p[i] = famJData->p[famJData->nParticles - 1];
                                /*
                                 * reset iterator and particle counter
                                 */
                                i--;
                                famJData->nParticles--;
                                quadData->nParticles++;
                        }
                }
#ifdef DEBUG
                printf("\n");
                printf("[pchase %i replace] replace done - resulting in quad(0x%08X,0x%08X)[%i]\n",
                       p4est->mpirank, p->x, p->y, quadData->nParticles);
#endif
        }
}

void
pchase_world_destroy(pchase_world_t * W)
{
        /* destroy allocated list */
        int                 i;
        sc_list_t          *tmp;
        sc_list_destroy(W->particle_push_list);
#ifdef DEBUG
        printf("[pchase %i world_destroy] destroyed particle_push_list\n", W->p4est->mpirank);
#endif
        for (i = 0; i < W->particles_to->elem_count; i++) {
                tmp = *((sc_list_t **) sc_array_index_int(W->particles_to, i));
                /* free every particle in tmp list */
                pchase_particle_t  *tmpParticle;
                while (tmp->first != NULL) {
                        tmpParticle = sc_list_pop(tmp);
                        P4EST_FREE(tmpParticle);
                }
                sc_list_destroy(tmp);
#ifdef DEBUG
                printf("[pchase %i world_destroy] destroyed particles_to[%i] list\n", W->p4est->mpirank, i);
#endif
        }
        sc_array_destroy(W->particles_to);
#ifdef DEBUG
        printf("[pchase %i world_destroy] destroyed particles_to_proc array\n", W->p4est->mpirank);
#endif
}
