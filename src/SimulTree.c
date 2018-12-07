#include "SimulTree.h"
#include <math.h>
#include "Utils.h"


#define INC_FOSSIL_ITEM 50
#define INFTY 9.9E99;

/* return  a random type of event wrt the rates*/
static char getType(gsl_rng *rgen, double birth, double death);

/* return  a random type of event wrt the rates*/
char getType(gsl_rng *rgen, double birth, double death) {
	double uni = gsl_rng_uniform(rgen);
	if(uni<birth/(birth+death))
		return 'b';
	else
		return 'd';
}



/*simulate a random tree with specified birth and death rates and fossil finds on this tree with rate "fossil"*/
TypeTree *simulTree(gsl_rng *rgen, double birth, double death, double maxTime) {
	TypeTree *tree;
	int *cur, ncur, i;
	double time = 0.;
	if((cur = (int*) malloc((MAX_CURRENT+1)*sizeof(int))) == NULL)
		return NULL;
	tree = newTree(INC_SIZE);
	tree->maxTime = maxTime;
	tree->maxTimeInt.inf = tree->maxTime;
	tree->maxTimeInt.sup = tree->maxTime;
	tree->minTime = 0.;
	tree->minTimeInt.inf = tree->minTime;
	tree->minTimeInt.sup = tree->minTime;
	tree->size = 1;
	tree->root = 0;
	tree->time[0] = INFTY;
	tree->node[0].child = -1;
	tree->node[0].sibling = -1;
	cur[0] = 0;
	ncur = 1;
	while(time<maxTime && ncur>0) {
		int type = getType(rgen, birth, death);
		int which = gsl_rng_uniform_int(rgen, ncur);
		double wait = gsl_ran_exponential(rgen, 1./(ncur*(birth+death)));
//		double wait = -log(gsl_ran_flat (rgen, 0., 1.))/(ncur*(birth+death+fossil));
		time += wait;
		if(time < maxTime) {
			switch(type) {
				case 'b':
					if(ncur > MAX_CURRENT) {
						printf("too much lineages generated during simulations (%d - max %d)\n", ncur, MAX_CURRENT);
						freeTree(tree);
						return NULL;
//						exit(EXIT_FAILURE);
					}
					tree->time[cur[which]] = time;
					if((tree->size+1)>=tree->sizeBuf) {
						tree->sizeBuf += INC_SIZE;
						tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
						tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
					}
					tree->node[cur[which]].child = tree->size;
					tree->time[tree->size] = INFTY;
					tree->node[tree->size].child = -1;
					tree->node[tree->size].sibling = tree->size+1;
					tree->time[tree->size+1] = INFTY;
					tree->node[tree->size+1].child = -1;
					tree->node[tree->size+1].sibling = -1;
					cur[which] = tree->size;
					cur[ncur] = tree->size+1;
					ncur++;
					tree->size += 2;
					break;
				case 'd':
					tree->time[cur[which]] = time;
					for(i=which+1; i<ncur; i++)
						cur[i-1] = cur[i];
					ncur--;
					break;
				default:
					break;
			}
		}
	}
	for(i=0; i<ncur; i++)
		tree->time[cur[i]] = maxTime;
	free((void*)cur);
	tree->sizeBuf = tree->size;
	if(tree->size) {
		tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
		tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
	} else {
		free((void*)tree->node);
		free((void*)tree->time);
		tree->node = NULL;
		tree->time = NULL;
	}
	return tree;
}
