#include "pchase_world.h"

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
        W->delta_t = 0.1;
        W->t_end = 1.0;
        W->n_particles = 0;
        W->step = 0;
        W->p4est = p4est;
        W->init_fn = init_fn;
        W->coarsen_fn = NULL;
        W->refine_fn = refine_fn;
        W->search_fn = search_fn;
        W->replace_fn = NULL;
        W->viter_fn = viter_fn;
        for (i = 0; i < DIM; i++)
                W->length[i] = 1.0;
        /* reset seed */
        srand(time(NULL));
#ifdef DEBUG
        /*
         * produce a random number and discard it. Otherwise rand() produces
         * almost everytime a number about 0.431..
         */
        rand();
#endif
        return W;
}

void
pchase_world_simulate(pchase_world_t * W)
{
        /* simulate until the end has come */
        while (W->t <= W->t_end) {
                pchase_world_update_x(W);
#ifdef DEBUG
                printf("Actual time on earth: %f\n", W->t);
#endif
                W->t += W->delta_t;
                W->step++;
        }
}

void
pchase_world_update_x(pchase_world_t * W)
{
        /* evaluate potential */
        printf("EVALUATE VELOCITY FIELD - NOT IMPLEMENTED YET");
}


pchase_particle_t  *
pchase_world_random_particle(pchase_world_t * W)
{
        int                 i;
        pchase_particle_t  *p = malloc((sizeof(pchase_particle_t)));
        /* pchase_particle_t  *p = P4EST_ALLOC(pchase_particle_t, 1); */

        for (i = 0; i < DIM; i++)
                p->x[i] = W->length[i] * rand() / (RAND_MAX + 1.);
#ifdef DEBUG
        printf("[pchase randPart] ");
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
        int                 i;

        /* get place for one point and take care of it via miniQuad */
        point = sc_array_new_size(sizeof(p4est_quadrant_t), 1);
        miniQuad = (p4est_quadrant_t *) sc_array_index(point, 0);

        /* create mini quadrant that is enclosing the given particle p */
        pchase_translate_particle_to_p4est(W, p, miniQuad);

#ifdef DEBUG
        printf("[pchase insertPart] Particle(x,y)=(%lf,%lf)\n", p->x[0], p->x[1]);
        printf("[pchase insertPart] in mini quad(%lld,%lld) at %lld\n", miniQuad->x, miniQuad->y, miniQuad);
#endif

        /*
         * check if the quadrant holding our to be inserted particle lies on
         * this processor
         */
        if (W->p4est->mpirank == 0) {
                /*
                 * find most deepest quadrant which encloses the mini quad in
                 * point and save its data in point.piggy3
                 */
                p4est_search(W->p4est, W->search_fn, point);

                /* extract enclosing quad from miniQuad.piggy3 */
                p4est_tree_t       *enclQuadTree = p4est_tree_array_index(W->p4est->trees, miniQuad->p.piggy3.which_tree);
                p4est_quadrant_t   *enclQuad = p4est_quadrant_array_index(&enclQuadTree->quadrants, miniQuad->p.piggy3.local_num);
                pchase_quadrant_data_t *enclQuadData = enclQuad->p.user_data;
#ifdef DEBUG
                /* print saved quad once again */
                printf("[pchase insertPart] ROOT is doing the job\n ");
                printf("[pchase insertPart] quad in user_pointer: %lld\n", enclQuad);
#endif

                /*
                 * TODO: - check if there are already 5 particles inside
                 * quads particle array and flag quad to refine - free all
                 * unneeded data - initialize quadData->nParticles in init_fn
                 */

#ifdef DEBUG
                printf("[pchase insertPart] #Particles in Quad: %d \n", enclQuadData->nParticles);
                printf("[pchase insertPart] BEFORE REAL INSERTION\n");
                /* print out all quadrants */
                p4est_iterate(W->p4est, NULL, NULL, W->viter_fn, NULL, NULL);
#endif
                /* reserve some memory for the particle struct */
                enclQuadData->p[enclQuadData->nParticles] = *P4EST_ALLOC(pchase_particle_t, 1);
                /* insert particle data into quad */
                enclQuadData->p[enclQuadData->nParticles].ID = p->ID;
                for (i = 0; i < DIM; ++i)
                        enclQuadData->p[enclQuadData->nParticles].x[i] = p->x[i];
                enclQuadData->nParticles++;

        } else
                /* send particle to its belonging proc */
                printf("[pchase insertPart] Not yet implemented");

#ifdef DEBUG
        printf("[pchase insertPart] REAL INSERTION DONE\n");
#endif
}

void
pchase_translate_particle_to_p4est(pchase_world_t * W, const pchase_particle_t * p, p4est_quadrant_t * q)
{
        double              quadrant_length = (double)(1 << P4EST_QMAXLEVEL);

#ifdef DEBUG
        printf("[pchase translatePart] Step 1: p->x = %lf\n", p->x[0]);
        printf("[pchase translatePart] Step 2: p->x/W->length.x = %lf\n", p->x[0] / W->length[0]);
        printf("[pchase translatePart] Step 3: times quad_len = %lf with quad_len = %lf and root_len = %lld\n", p->x[0] / W->length[0] * quadrant_length, quadrant_length, P4EST_ROOT_LEN);
        printf("[pchase translatePart] Step 4: truncate = %d\n", (p4est_qcoord_t) (p->x[0] / W->length[0] * quadrant_length));
        printf("[pchase translatePart] Step 5: times 2 gives q->x = %lld\n", (p4est_qcoord_t) (p->x[0] / W->length[0] * quadrant_length) << 1);
#endif

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
        p4est_quadrant_t   *q;
        int                 quadrant_length;

#ifdef DEBUG
        printf("[pchase search] HELLOOO I AM A CALLBACK in x:%lld y:%lld with tree:%lld and local_num:%lld\n", quadrant->x, quadrant->y, which_tree, quadrant->p.piggy3.local_num);
#endif

        quadrant_length = P4EST_QUADRANT_LEN(quadrant->level);
        q = (p4est_quadrant_t *) point;

        /* mini quad must lie entirely in quadrant */
        if (q->x >= quadrant->x && q->x <= quadrant->x + quadrant_length &&
            q->y >= quadrant->y && q->y <= quadrant->y + quadrant_length) {
                /* only save the quad in case it's a leaf */
                if (is_leaf) {
                        /* save current quad via piggy3 */
                        q->p.piggy3.local_num = quadrant->p.piggy3.local_num;
                        q->p.piggy3.which_tree = which_tree;
                }
#ifdef DEBUG
                printf("[pchase search] YES YES YES - we found a quadrant whos child holds " \
                       "our mini quad at linear positon: %lld in level: %lld is_leaf: %d\n",
                       p4est_quadrant_linear_id(quadrant, quadrant->level), quadrant->level, is_leaf);
                printf("[pchase search] QuadLen: %lld, pos(x,y)=(%lld,%lld)\n", P4EST_QUADRANT_LEN(quadrant->level), quadrant->x, quadrant->y);
                printf("[pchase search] Found this quad pointer: %lld\n", quadrant);
#endif

                return 1;
        } else
                return 0;
}

static void
init_fn(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
        printf("[pchase init_fn] BUM, INIT for Quad(x,y): (%lld,%lld) at memory: %lld\n", quadrant->x, quadrant->y, quadrant);

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
        pchase_quadrant_data_t *qData = (pchase_quadrant_data_t *) info->quad->p.user_data;
        for (i = 0; i < qData->nParticles; i++) {
                printf("[pchase destroy_fn] FREEEE, free like the wind\n");
                P4EST_FREE(&qData->p[i]);
        }
        return;
}
static void
viter_fn(p4est_iter_volume_info_t * info, void *user_data)
{
        pchase_quadrant_data_t *quadData = (pchase_quadrant_data_t *) info->quad->p.user_data;
        printf("[pchase main iterate] quad(%lld,%lld) with %i particles \n", info->quad->x, info->quad->y, quadData->nParticles);
        return;
}
