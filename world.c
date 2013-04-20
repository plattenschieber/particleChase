#include "world.h"
#include "world.h"
#include <stdio.h>
#include <stdlib.h>

void simulate(LIST *L) {
	while(world_t < world_t_end) {
		update_x();
		world_t+=world_delta_t;
	}
}
void update_x() {
	/* do some stuff */
}

