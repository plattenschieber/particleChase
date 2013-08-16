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
        p->ID = W->n_particles;
#endif
        W->n_particles++;

        return p;
}

p4est_quadrant_t   *
pchase_world_insert_particle(pchase_world_t * W, pchase_particle_t * p)
{
        /* get mini quadrant in which the particle lies */
        p4est_quadrant_t   *q = pchase_translate_particle_to_p4est(W, p);

        /*
         * check if the quadrant holding our to be inserted particle lies on
         * this proc
         */
        p4est_comm_find_owner(W->p4est, W->p4est->first_local_tree, q, -1);
        scanf("%i", NULL);
        if (1 == W->p4est->mpirank) {
                printf("Not yet implemented");
        }
        /* send particle to belonging */
        else
                printf("Not yet implemented");

        return q;
}

p4est_quadrant_t   *
pchase_translate_particle_to_p4est(pchase_world_t * W, pchase_particle_t * p)
{
        p4est_quadrant_t   *q;
        q = P4EST_ALLOC(p4est_quadrant_t, 1);
        q->level = P4EST_MAXLEVEL;
        /* first convert particles' coord to [0,1) then place it into p4est */
        q->x = (int)floor((p->x[0] / W->length[0]) * P4EST_ROOT_LEN);
        q->y = (int)floor((p->x[1] / W->length[1]) * P4EST_ROOT_LEN);
#if DIM == 3
        q->z = (int)floor((p->x[2] / W->length[2]) * P4EST_ROOT_LEN);
#endif
        return q;
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
