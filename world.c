#include "world.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

WORLD *world_init(FILE *parameter, FILE *particles) {
	/* get some place for the world */
	WORLD *W;
	if (!(W=malloc(sizeof(WORLD)))) {
		fprintf(stderr, "Could not allocate memory for the World\n");
		exit(EXIT_FAILURE); 
	}
	/* get some place for its cells */
	W->cells = list_init();
	/* read in and set all parameters */
	W->n_particles = 0;
	W->step = 0;
	world_read_parameter(W, parameter);
	/* world_read_particles(W, particles); */
	/* reset seed */
	srand(time(NULL));
	/* fill in some particles randomly instead */
	world_randomfill(W);
	return W;
}

void world_simulate(WORLD *W) {
	/* simulate until the end has come */
	while (W->t <= W->t_end) {
		world_update_x(W);
#ifdef DEBUG
		printf("Actual time on earth: %f\n",W->t);
#endif
		print_particles(W->particles);
		W->t+=W->delta_t;
		W->step++;
	}
}

void world_update_x(WORLD *W) {
	NODE *it;
	PARTICLE *p;
	double x,y;
 	for (it = W->particles->first; it != NULL; it = it->next) {
		p = it->data;
		x = -p->x[1];
		y = p->x[0];

		double norm = sqrt(x*x + y*y);
		p->x[0] += W->delta_t * x/norm;
		p->x[1] += W->delta_t * y/norm;
	}
}

void world_update_v(WORLD *W) {
	return;
}

void world_read_particles(WORLD *W, FILE *particles) {
	int ID, i;
	PARTICLE *tmp;
	while (fscanf(particles,"%i ",&ID) != EOF) { 
		/* get new memory for a temporary 'PARTICLE' */
		tmp = malloc(sizeof(PARTICLE)); 
		tmp->ID = ID;
#ifdef DEBUG
		printf("\n%i\t",tmp->ID);
#endif
		/* read in the coordinates */ 
		for (i=0; i<DIM; i++) {
			if(fscanf(particles,"%lf", &tmp->x[i])!=1) {
				fprintf(stderr, "Something is wrong with your particle file - Abort\n");
				exit(EXIT_FAILURE);
			}	
#ifdef DEBUG
			printf("%f\t",tmp->x[i]);
#endif
		}

		list_insert(W->particles, NULL, tmp);
		W->n_particles++; 
#ifdef DEBUG
		printf("particle inserted\n");
#endif
	}
#ifdef DEBUG
	printf("Read %i Particle(s)\n", W->n_particles);
#endif
	
}
void world_read_parameter(WORLD *W, FILE *parameter) {
	char *option = malloc(sizeof(char)*20);
	int i,line=0;
#ifdef DEBUG
	printf("read parameter\n");
#endif
	while(fscanf(parameter,"%s", option) == 1)
	{	
		/* update current linenumber */ 
		line++;
#ifdef DEBUG
		printf("%s\n",option);
#endif
		if (strcmp(option,"delta_t") == 0)
			fscanf(parameter, "%lf", &W->delta_t);
		else if (strcmp(option,"t_end") == 0)
			fscanf(parameter, "%lf", &W->t_end);
		else if (strcmp(option,"t") == 0)
			fscanf(parameter, "%lf", &W->t);
		else if (strcmp(option,"length") == 0)
			for (i=0; i<DIM; i++)
				fscanf(parameter, "%lf", &W->length[i]);
		else if (strcmp(option,"step") == 0)
			fscanf(parameter, "%i", &W->step);
		else if (strcmp(option,"n_particles") == 0)
			fscanf(parameter, "%i", &W->n_particles);
		else {
			printf("Could not read option in line %i of given parameter file. Abort.\n", line);
			exit(EXIT_FAILURE);
		}
	}
}

void world_randomfill(WORLD *W) {
	unsigned i;
	for (i=0; i<MAX_PARTICLES; i++) 
		world_random_particle(W);

}

PARTICLE * random_particle() {
	int i;
	PARTICLE *p = malloc(sizeof(PARTICLE));
	for (i=0; i<DIM; i++)
		p->x[i] = rand()/(RAND_MAX + 1.);
	p->ID = PARTICLE_COUNT;
	PARTICLE_COUNT++; 
	return p;
}
PARTICLE * world_random_particle(WORLD *W) {
	int i;
	PARTICLE *p = malloc(sizeof(PARTICLE));
	for (i=0; i<DIM; i++)
		p->x[i] = W->length[i] * rand()/(RAND_MAX + 1.);
	p->ID = W->n_particles;
	list_insert(W->particles, NULL, p);
	W->n_particles++; 
	return p;
}

void world_print_particlesXYZ(WORLD *W) {
	int i;
	NODE *it;
	PARTICLE *p;
	printf("%i\n",W->n_particles);
	//for (it = L->first; it; it = it->next) {
		/* cast data to particle */
		p = it->data;
		printf("%i\t",p->ID);
		for (i = 0; i < DIM; i++)
			printf("%f\t",p->x[i]);
		printf("\n");
	//}
	/*
    int nParticles = 0;
    // write size and actual time of our world W
    //coordinates << W.nParticles << std::endl << "Time: " << W.t << std::endl;
    for (std::vector<Cell>::const_iterator i = W.particles.begin(); i < W.particles.end(); i++)
        for (std::list<Particle>::const_iterator j = i->particles.begin(); j != i->particles.end(); j++)
            nParticles++;
    coordinates << nParticles << std::endl << "Time: " << W.t << std::endl;
    // get out every particle to satisfy the xyz format
    for (std::vector<Cell>::const_iterator i = W.particles.begin(); i < W.particles.end(); i++)
        for (std::list<Particle>::const_iterator j = i->particles.begin(); j != i->particles.end(); j++)
        {
            // each particle should be an H-atom. At least now...
            coordinates << "H\t";
            // particle j is located in a DIM-dimensional space
            for (unsigned int d=0; d<DIM; d++)
                // get it out, seperated with tabulars
                coordinates << j->x[d] << "\t";
            // newline at end of each particle
            coordinates << std::endl;
        }
*/

}
