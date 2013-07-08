#include <stdio.h>
#include <stdlib.h>
#include "particle.h"
#include "world.h"

int main() {
	int i;
	PARTICLE p1;
   	PARTICLE p2;
	LIST *L = list_init();		
	NODE *it;

	/* insert some particles */
	p1.x[0] = 1.3; p1.x[1] = 2.3;
	p1.ID = 2;

	p2.x[0] = 1.5; p2.x[1] = 4.8;
	p2.ID = 1;
	list_insert(L, L->first, &p2);

	/* print particles before and after the simulation */ 
	print_particles(L);
	simulate(L);
	print_particles(L);

	return 0;
}
