#include <stdio.h>
#include <stdlib.h>
#include "particle.h"

int main() {
	int i;
	PARTICLE particle;
	LIST *L;
	particle.x[0] = 1.3; particle.x[1] = 2.3;
	particle.ID = 2;
	L->first = list_create(L, &particle);

	particle.x[0] = 1.5; particle.x[1] = 4.8;
	particle.ID = 1;
	list_insert(L, L->first, &particle);

	}
	return 0;
}
