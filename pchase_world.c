#include "pchase_world.h"

static FILE * pchase_output;

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
        W->delta_t = 0.00001;
        W->t_end = 9.00;
        W->n_particles = 0;
        W->step = 0;
        W->p4est = p4est;
        W->init_fn = init_fn;
        W->coarsen_fn = NULL;
        W->refine_fn = refine_fn;
        W->search_fn = search_fn;
        W->replace_fn = NULL;
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
        /* simulate until the end has come */
        while (W->t <= W->t_end) {
#ifdef PRINTXYZ
                if(W->step % 1000 == 0) {
                        fprintf(pchase_output,"%i\n",W->n_particles);    
                        fprintf(pchase_output,"Time: %f\n", W->t);
                }
#endif
                p4est_iterate(W->p4est, NULL, W, W->update_x_fn, NULL, NULL);
                W->t += W->delta_t;
                W->step++;
        }
        printf("simulation over\n");
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
        printf("[pchase %i insertPart] Translated Particle(%lf,%lf) to miniQuad (0x%08X,0x%08X)\n",
             W->p4est->mpirank, p->x[0], p->x[1], miniQuad->x, miniQuad->y);
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

                /* insert particle data into quad and update particle counter */
                enclQuadData->p[enclQuadData->nParticles] = p;
                enclQuadData->nParticles++;
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
                printf("[pchase %i insertPart] freed particle\n", W->p4est->mpirank);
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
        printf("[pchase %i init_fn] quad(0x%08X,0x%08X) in tree: %lld\n", p4est->mpirank, quadrant->x, quadrant->y, which_tree);
        ((pchase_quadrant_data_t *) quadrant->p.user_data)->nParticles = 0;
}

static int
refine_fn(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
        if (0)
                return 1;
        else
                return 0;
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
#elif
                printf("%lf\t%lf\n", quadData->p[i]->x[0], quadData->p[i]->x[1]);
#endif
        
}
static void
update_x_fn(p4est_iter_volume_info_t * info, void *user_data)
{
        double              x, y, norm;
        int                 i;
        pchase_quadrant_data_t *quadData = (pchase_quadrant_data_t *) info->quad->p.user_data;
        pchase_world_t * W = (pchase_world_t *) user_data;


        for (i = 0; i < quadData->nParticles; i++) {
                x = -quadData->p[i]->x[1] + 0.5;
                y = quadData->p[i]->x[0] - 0.5;

                norm = sqrt(x * x + y * y);
                quadData->p[i]->x[0] += W->delta_t * x / norm;
                quadData->p[i]->x[1] += W->delta_t * y / norm;
                if(W->step % 100 == 0)
#ifdef PRINTXYZ
                        fprintf(pchase_output,"H\t%lf\t%lf\t0\n", quadData->p[i]->x[0], quadData->p[i]->x[1]);
#elif defined(PRINTGNUPLOT)
                        fprintf(pchase_output,"%lf\t%lf\n", quadData->p[i]->x[0], quadData->p[i]->x[1]);
#endif
                /* move particle if it has left the quad */ 
                if (!pchase_particle_lies_in_quad(W,quadData->p[i],info->quad)) {
                        pchase_world_insert_particle(W, quadData->p[i]);

                        /* if it's not the last particle in the array */
                        if (i!=quadData->nParticles-1) {
                                /* move the last particle to i'th place to prevent a hole */
                                quadData->p[i] = quadData->p[quadData->nParticles-1];
                                quadData->p[quadData->nParticles-1] = NULL;
                                /* set iterator accordingly */
                                i--;
                        }
                        /* update particle counter in quad array */
                        quadData->nParticles--;
                }
        }
}

int
pchase_particle_lies_in_quad(pchase_world_t * W, const pchase_particle_t * p, p4est_quadrant_t * q)
{
        p4est_qcoord_t quadrant_length,root_len;
        quadrant_length = P4EST_QUADRANT_LEN(q->level);
        root_len = P4EST_ROOT_LEN;

        if (p->x[0]*root_len < q->x || p->x[0]*root_len >= q->x + quadrant_length ||
            p->x[1]*root_len < q->y || p->x[1]*root_len >= q->y + quadrant_length) {
#ifdef DEBUG
                printf("particles at p4est_coord(0x%08X, 0x%08X) left Quad(%09lld,%lld)\n",
                                (int) (p->x[0]*root_len), (int) (p->x[1]*root_len), q->x, q->y);
#endif
                return 0;
        }
        else
                return 1;
}
