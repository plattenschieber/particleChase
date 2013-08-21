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
        W->p4est = p4est;
        W->init_fn = init_fn;
        W->coarsen_fn = coarsen_fn;
        W->refine_fn = refine_fn;
        W->search_fn = search_fn;
        W->replace_fn = replace_fn;
        W->viter_fn = viter_fn;
        W->destroy_fn = destroy_fn;
        W->print_fn = print_fn;
        W->update_x_fn = update_x_fn;
#ifdef PRINTXYZ
        pchase_output = fopen("pchase_particles.xyz", "w");
#elif defined(PRINTGNUPLOT)
        pchase_output = fopen("pchase_particles.plot", "w");
#endif
        for (i = 0; i < DIM; i++)
                W->length[i] = 1.0;
        /* reset seed */
        srand(time(NULL));
        return W;
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
        /* simulate until the end has come */
        while (W->t <= W->t_end) {
#ifdef PRINTXYZ
                if (W->step % 1000 == 0) {
                        fprintf(pchase_output, "%i\n", W->n_particles);
                        fprintf(pchase_output, "Time: %f\n", W->t);
                }
#endif
                /* update the position of all particles on all quads */
                p4est_iterate(W->p4est, NULL, W, W->update_x_fn, NULL, NULL);

                if (W->step % 1 == 0) {
                        /* refine every quad containing more than 5 particles */
                        p4est_refine_ext(W->p4est, 0, -1, W->refine_fn, W->init_fn, W->replace_fn);
                        p4est_coarsen_ext(W->p4est, 0, W->coarsen_fn, W->init_fn, W->replace_fn);
                        p4est_balance_ext(W->p4est, P4EST_CONNECT_FULL, W->init_fn, W->replace_fn);
                        p4est_partition_ext(W->p4est, 1, NULL);

                        /*
                         * convert current step to filename and write VTK
                         * Data entry
                         */
                        sprintf(fileNumber, "%03d", W->step);
                        strcat(VTKData, VTKData1);
                        strcat(VTKData, fileNumber);
                        strcat(VTKData, VTKData2);
                        strcat(fileName, fileName1);
                        strcat(fileName, fileNumber);
                        strcat(VTKData, fileName);
                        strcat(VTKData, VTKData3);
                        fprintf(vtk_timeseries, VTKData);
                        p4est_vtk_write_file(W->p4est, NULL, fileName);

                        i++;
                        /* reset filenames */
                        VTKData[0] = '\0';
                        fileName[0] = '\0';

                }
                W->t += W->delta_t;
                W->step++;
        }
        printf("simulation over\n");
        fprintf(vtk_timeseries, "      </Collection>\n");
        fprintf(vtk_timeseries, "</VTKFile>\n");
        fclose(vtk_timeseries);
        fclose(pchase_output);
}

pchase_particle_t  *
pchase_world_random_particle(pchase_world_t * W)
{
        int                 i;
        pchase_particle_t  *p = P4EST_ALLOC(pchase_particle_t, 1);

        for (i = 0; i < DIM; i++)
                p->x[i] = W->length[i] * rand() / (RAND_MAX + 1.);
#ifdef DEBUG
        printf("[pchase %i randPart] ", W->p4est->mpirank);
        p->ID = W->n_particles;
        for (i = 0; i < DIM; i++)
                printf("p->x[%d]=%lf,\t", i, p->x[i]);
        printf("\n");
#endif
        W->n_particles++;
        return p;
}

void
pchase_world_insert_particle(pchase_world_t * W, pchase_particle_t * p)
{
        p4est_quadrant_t   *miniQuad;
        sc_array_t         *point;
        int                 owner;

        /* get place for one point and take care of it via miniQuad */
        point = sc_array_new_size(sizeof(p4est_quadrant_t), 1);
        miniQuad = (p4est_quadrant_t *) sc_array_index(point, 0);

        /* create mini quadrant that is enclosing the given particle p */
        pchase_translate_particle_to_p4est(W, p, miniQuad);

#ifdef DEBUG
        printf("[pchase %i insertPart] Translated Particle[%i](%lf,%lf) to miniQuad (0x%08X,0x%08X)\n",
               W->p4est->mpirank, p->ID, p->x[0], p->x[1], miniQuad->x, miniQuad->y);
#endif
        /*
         * find most deepest quadrant which encloses the mini quad in point
         * and save its data in point.piggy3
         */
        p4est_search(W->p4est, W->search_fn, point);



        /* if the miniQuad is not flagged, it's lying on this proc */
        if (miniQuad->p.piggy3.local_num != -1) {
                /* extract enclosing quad from miniQuad.piggy3 */
                p4est_tree_t       *enclQuadTree = p4est_tree_array_index(W->p4est->trees, miniQuad->p.piggy3.which_tree);
                p4est_quadrant_t   *enclQuad = p4est_quadrant_array_index(&enclQuadTree->quadrants, miniQuad->p.piggy3.local_num);
                pchase_quadrant_data_t *enclQuadData = enclQuad->p.user_data;

                if (enclQuadData->nParticles < 25) {
                        /*
                         * insert particle data into quad and update particle
                         * counter
                         */
                        enclQuadData->p[enclQuadData->nParticles] = p;
                        enclQuadData->nParticles++;
                }
                /* too many particles in quad */
                else {
                        printf("[pchase %i insertPart] Too many (%i) particles ", W->p4est->mpirank, enclQuadData->nParticles);
                        if (enclQuad->level < P4EST_QMAXLEVEL) {
                                printf("- refining enclQuad and inserting particle afterwards\n");
                                p4est_refine_ext(W->p4est, 0, -1, W->refine_fn, W->init_fn, W->replace_fn);
                                pchase_world_insert_particle(W, p);
                        } else {
                                printf("- enclQuad is not refinable - we have to dissmiss this particle\n");
                                P4EST_FREE(p);
                        }
                }
#ifdef DEBUG
                /* print number of particles in quad */
                printf("[pchase %i insertPart] #Particles in enclQuad: %d \n", W->p4est->mpirank, enclQuadData->nParticles);
#endif
        }
        /* particle lies on another proc */
        else {
                /* send particle to its belonging proc */
#ifdef DEBUG
                printf("[pchase %i insertPart] Sending Particle not implemented yet\n", W->p4est->mpirank);
#endif

                /* use find_owner to pigeon-hole particle into the right proc */
                owner = p4est_comm_find_owner(W->p4est, miniQuad->p.piggy3.which_tree, miniQuad, W->p4est->mpirank);
#ifdef DEBUG
                printf("[pchase %i insertPart] particle[%i] sent to proc %i\n", W->p4est->mpirank, p->ID, owner);
#endif
                /* and free particle on this proc */
                P4EST_FREE(p);
#ifdef DEBUG
                printf("[pchase %i insertPart] freed particle\n", W->p4est->mpirank);
#endif
        }
        sc_array_destroy(point);
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
        if (miniQuad->x >= quadrant->x && miniQuad->x <= quadrant->x + quadrant_length &&
            miniQuad->y >= quadrant->y && miniQuad->y <= quadrant->y + quadrant_length) {
                /* remember current tree (needed to identify proc) */
                miniQuad->p.piggy3.which_tree = which_tree;
                /*
                 * if it's a leaf, it holds it's position in the tree and can
                 * only exist on this proc.
                 */
                if (is_leaf)
                        miniQuad->p.piggy3.local_num = quadrant->p.piggy3.local_num;
                /* flag the quad in the other case */
                else
                        miniQuad->p.piggy3.local_num = -1;
#ifdef DEBUG
                printf("[pchase %i search] found quad (holding miniQuad) " \
                       "with local_num %d in level: %d is_leaf: %d\n", p4est->mpirank,
                    miniQuad->p.piggy3.local_num, quadrant->level, is_leaf);
#endif

                return 1;
        } else
                return 0;
}

static void
init_fn(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
#ifdef DEBUG
        printf("[pchase %i init_fn] quad(0x%08X,0x%08X) in tree: %lld\n", p4est->mpirank, quadrant->x, quadrant->y, which_tree);
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
destroy_fn(p4est_iter_volume_info_t * info, void *Data)
{
        int                 i;
        pchase_particle_t  *p;
        pchase_quadrant_data_t *qData = (pchase_quadrant_data_t *) info->quad->p.user_data;
        for (i = 0; i < qData->nParticles; i++) {
                p = qData->p[i];
                P4EST_FREE(p);
#ifdef DEBUG
                printf("[pchase %i destroy_fn] freed particle[%i](%lf,%lf)\n", info->p4est->mpirank, p->ID, p->x[0], p->x[1]);
#endif
        }
        return;
}
static void
viter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
        pchase_quadrant_data_t *quadData = (pchase_quadrant_data_t *) info->quad->p.user_data;
        printf("[pchase %i main iterate] quad(0x%08X,0x%08X) has %i particles \n", info->p4est->mpirank, info->quad->x, info->quad->y, quadData->nParticles);
        return;
}
static void
print_fn(p4est_iter_volume_info_t * info, void *user_data)
{
        int                 i;
        pchase_quadrant_data_t *quadData = (pchase_quadrant_data_t *) info->quad->p.user_data;

        for (i = 0; i < quadData->nParticles; i++)
#ifdef DEBUG
                printf("%i\t%lf\t%lf\n", quadData->p[i]->ID, quadData->p[i]->x[0], quadData->p[i]->x[1]);
#else
                printf("%lf\t%lf\n", quadData->p[i]->x[0], quadData->p[i]->x[1]);
#endif

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
}

static void
update_x_fn(p4est_iter_volume_info_t * info, void *user_data)
{
        int                 i;
        pchase_quadrant_data_t *quadData = (pchase_quadrant_data_t *) info->quad->p.user_data;
        pchase_world_t     *W = (pchase_world_t *) user_data;

        for (i = 0; i < quadData->nParticles; i++) {
                /* update particles' velocity */
                pchase_world_velocity(W, quadData->p[i]);

#ifdef DEBUG
                printf("[pchase %i updateX] particle[%i](%lf,%lf) in quad(%lld) with %i particles\n",
                       info->p4est->mpirank, quadData->p[i]->ID, quadData->p[i]->x[0], quadData->p[i]->x[1], info->quadid, quadData->nParticles);
#endif
#ifdef PRINTXYZ
                fprintf(pchase_output, "H\t%lf\t%lf\t0\n", quadData->p[i]->x[0], quadData->p[i]->x[1]);
#elif defined(PRINTGNUPLOT)
                fprintf(pchase_output, "%lf\t%lf\n", quadData->p[i]->x[0], quadData->p[i]->x[1]);
#endif
                /* move particle if it has left the quad */
                if (!pchase_particle_lies_in_quad(quadData->p[i], info->quad)) {
#ifdef DEBUG
                        printf("[pchase %i updateX] particle[%i](%lf,%lf) left quad(0x%08X,0x%08X)\n",
                               info->p4est->mpirank, quadData->p[i]->ID, quadData->p[i]->x[0], quadData->p[i]->x[1], info->quad->x, info->quad->y);
#endif
                        pchase_world_insert_particle(W, quadData->p[i]);

                        /* if it's not the last particle in the array */
                        if (i != quadData->nParticles - 1) {
                                /*
                                 * move the last particle to i'th place to
                                 * prevent a hole
                                 */
                                quadData->p[i] = quadData->p[quadData->nParticles - 1];
#ifdef DEBUG
                                quadData->p[quadData->nParticles - 1] = NULL;
#endif
                                /* set iterator accordingly */
                                i--;
                        }
                        /* update particle counter in quad array */
                        quadData->nParticles--;
                }
        }
}

int
pchase_particle_lies_in_quad(const pchase_particle_t * p, p4est_quadrant_t * q)
{
        p4est_qcoord_t      quadrant_length, root_len;
        quadrant_length = P4EST_QUADRANT_LEN(q->level);
        root_len = P4EST_ROOT_LEN;

        if (p->x[0] * root_len < q->x || p->x[0] * root_len >= q->x + quadrant_length ||
            p->x[1] * root_len < q->y || p->x[1] * root_len >= q->y + quadrant_length) {
#ifdef DEBUG
                printf("particle[%i] at p4est_coord(%09lld,%09lld) left Quad(0x%08X, 0x%08X)\n",
                       p->ID, (int)(p->x[0] * root_len), (int)(p->x[1] * root_len), q->x, q->y);
#endif
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
                printf("[pchase %i replace] REPLACING QUAD BY %i CHILDREN\n", p4est->mpirank, P4EST_CHILDREN);
#endif
                /* set readable names */
                p = outgoing[0];
                quadData = (pchase_quadrant_data_t *) p->p.user_data;
                fam = incoming;
                for (i = 0; i < quadData->nParticles; i++)
                        for (j = 0; j < P4EST_CHILDREN; j++) {
                                if (pchase_particle_lies_in_quad(quadData->p[i], fam[j])) {
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
#ifdef DEBUG
                                        quadData->p[quadData->nParticles - 1] = NULL;
#endif
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
                printf("[pchase %i replace] REPLACING %i CHILDREN BY THEIR PARENT\n", p4est->mpirank, P4EST_CHILDREN);
#endif
                /* set readable names */
                p = incoming[0];
                quadData = (pchase_quadrant_data_t *) p->p.user_data;
                fam = outgoing;
                for (j = 0; j < P4EST_CHILDREN; j++) {
                        famJData = (pchase_quadrant_data_t *) fam[j]->p.user_data;
                        for (i = 0; i < famJData->nParticles; i++) {
                                /* move particle to parent */
                                quadData->p[quadData->nParticles] = famJData->p[i];
                                /*
                                 * move last particle in array to originated
                                 * hole
                                 */
                                famJData->p[i] = famJData->p[famJData->nParticles - 1];
#ifdef DEBUG
                                famJData->p[famJData->nParticles - 1] = NULL;
#endif
                                /* reset iterator and particle counter */
                                i--;
                                famJData->nParticles--;
                                quadData->nParticles++;
                        }
                }
        }
}
