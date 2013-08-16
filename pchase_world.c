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
        pchase_particle_t  *p = P4EST_ALLOC(pchase_particle_t, 1);

        for (i = 0; i < DIM; i++)
                p->x[i] = W->length[i] * rand() / (RAND_MAX + 1.);
#ifdef DEBUG
        printf("Random Particle Generator: ");
        p->ID = W->n_particles;
        for (i = 0; i < DIM; i++)
                printf("p->x[%d]=%lf,\t", i, p->x[i]);
        printf("\n");
#endif
        W->n_particles++;
        return p;
}

p4est_quadrant_t   *
pchase_world_insert_particle(pchase_world_t * W, pchase_particle_t * p)
{
        p4est_quadrant_t   *q;
        sc_array_t         *point;

        /* get place for one point and take care of it via q */ 
        point = sc_array_new_size(sizeof(p4est_quadrant_t), 1); 
        q = (p4est_quadrant_t *) sc_array_index(point, 0);

        /* actually create mini quadrant that is enclosing the given partile p */
        pchase_translate_particle_to_p4est(W, p, q);

#ifdef DEBUG
        printf("Insert Particle: particle(x,y)=(%lf,%lf)\n", p->x[0], p->x[1]);
        printf("quadrant(x,y)=(%d,%d)\n", q->x, q->y);
#endif

        /*
         * check if the quadrant holding our to be inserted particle lies on
         * this processor
         * 
         * p4est_comm_find_owner(W->p4est, W->p4est->first_local_tree, q, * -1);
         */
        if (W->p4est->mpirank == 0) {
                p4est_search(W->p4est, W->search_fn, point);

        }
        /* send particle to belonging */
        else
                printf("Not yet implemented");

        return q;
}

void
pchase_translate_particle_to_p4est(pchase_world_t * W, const pchase_particle_t * p, p4est_quadrant_t * q)
{
        double              quadrant_length = (double)(1 << P4EST_QMAXLEVEL);

#ifdef DEBUG
        printf("Step 1: p->x = %lf\n", p->x[0]);
        printf("Step 2: p->x/W->length.x = %lf\n", p->x[0] / W->length[0]);
        printf("Step 3: times quad_len = %lf with quad_len = %lf and root_len = %lld\n", p->x[0] / W->length[0] * quadrant_length, quadrant_length, P4EST_ROOT_LEN);
        printf("Step 4: truncate = %d\n", (p4est_qcoord_t) (p->x[0] / W->length[0] * quadrant_length));
        printf("Step 5: q->x = %d\n", (p4est_qcoord_t) (p->x[0] / W->length[0] * quadrant_length) << 1);
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
        printf("HELLOOO I AM A CALLBACK in x:%i y%i\n", quadrant->x, quadrant->y);
#endif

        quadrant_length = P4EST_QUADRANT_LEN(quadrant->level);
        q = (p4est_quadrant_t *) point;

        if (q->x >= quadrant->x && q->x <= quadrant->x + quadrant_length &&
            q->y >= quadrant->y && q->y <= quadrant->y + quadrant_length &&
            /* q->p.user_int < quadrant->level) { */
            1) {
#ifdef DEBUG
                printf("YES YES YES - we found a quadrant whos child holds our mini quad\n");
#endif 

                /*
                 * replace current_quad with quadrant and save it's level for
                 * faster access in user_int
                 */
                q->p.user_int = quadrant->level;
                return 1;
        }
        else 
            return 0;
}

static void
init_fn(p4est_t * p4est, p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
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
