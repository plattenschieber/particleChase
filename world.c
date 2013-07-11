#include "world.h"
#include <stdio.h>
#include <stdlib.h>

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
	world_read_particles(W, particles);
	world_read_parameter(W, parameter);
	return W;
}

void world_simulate(WORLD *W) {
	/* simulate until the end has come */
	while (W->t <= W->t_end) {
		world_update_x(W);
#ifdef DEBUG
		printf("Actual time on earth: %f\n",W->t);
#endif
		W->t+=W->delta_t;
		W->step++;
	}
}

void world_update_x(WORLD *W) {
	NODE *it;
	PARTICLE *p;
	unsigned int i;
 	for (it = W->cells->first; it != NULL; it = it->next)
		for (i=0; i<DIM; i++) {
			/* cast data to particle */
			p = it->data;
			/* do some calculation (move more in y coordinate than x)*/
			p->x[i] += 0.112+i/2;
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
#ifdef DEBUG
			printf("%f\t",tmp->x[i]);
#endif
		}

		list_insert(W->cells, NULL, tmp);
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
			printf("Could not read option in line %i of given parameter file. Abort.", line);
			exit(EXIT_FAILURE);
		}
	}
}
