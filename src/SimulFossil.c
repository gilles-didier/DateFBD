#include <math.h>
#include "SimulFossil.h"
#include "Utils.h"


#define INC_FOSSIL_ITEM 50
#define INFTY 9.9E99;


static char getType(gsl_rng *rgen, double birth, double death, double fossil);


/* return  a random type of event wrt the rates*/
char getType(gsl_rng *rgen, double birth, double death, double fossil) {
	double uni = gsl_rng_uniform(rgen);
	if(uni<birth/(birth+death+fossil)) {
		return 'b';
	} else {
			if(uni<(birth+death)/(birth+death+fossil))
			return 'd';
		else
			return 'f';
	}
}


/*adds fossil finds with rate "fossil" on the tree "tree"*/
TypeFossilFeature *addFossils(gsl_rng *rgen, double fossil, TypeTree *tree) {
	int cur[MAX_CURRENT], ncur, i, n;
	double time = 0., maxTime;
	TypeFossilFeature *fos;
	fos = (TypeFossilFeature*) malloc(sizeof(TypeFossilFeature));
	fos->status = NULL;
	fos->sizeBuf = INC_FOSSIL_ITEM;
	fos->fossil = (int*) malloc(tree->size*sizeof(int));
	for(n=0; n<tree->size; n++)
		fos->fossil[n] = -1;
	fos->fossilList = (TypeFossilList*) malloc(fos->sizeBuf*sizeof(TypeFossilList));
	fos->size = 0;
	if(tree->time[tree->root] == 0.) {
		ncur = 2;
		cur[0] = tree->node[tree->root].child;
		cur[1] = tree->node[tree->node[tree->root].child].sibling;
	} else {
		ncur = 1;
		cur[0] = tree->root;
	}
	maxTime = 0.;
	for(i=0; i<tree->size; i++)
		if(tree->time[i]>maxTime)
			maxTime = tree->time[i];
	time = 0.;
	while(ncur>0 && time<maxTime) {
		double next;
		int lin = 0;
		next = tree->time[cur[0]];
		for(i=1; i<ncur; i++) {
			if(tree->time[cur[i]] < next) {
				next = tree->time[cur[i]];
				lin = i;
			}
		}
		do {
			time += gsl_ran_exponential(rgen, 1./(ncur*fossil));
			if(time <= next) {
				int which = gsl_rng_uniform_int(rgen, ncur);
				if(fos->size>=fos->sizeBuf) {
					fos->sizeBuf += INC_SIZE;
					fos->fossilList = (TypeFossilList*) realloc((void*)fos->fossilList, fos->sizeBuf*sizeof(TypeFossilList));
				}
				fos->fossilList[fos->size].time = time;
				fos->fossilList[fos->size].prec = fos->fossil[cur[which]];
				fos->fossil[cur[which]] = fos->size;
				fos->size++;
			}
		} while(time<=next);
		if(tree->node[cur[lin]].child>=0 && tree->node[tree->node[cur[lin]].child].sibling>=0) {
			if(ncur >= MAX_CURRENT) {
				printf("too much lineages generated during simulations");
				exit(EXIT_FAILURE);
			}
			cur[ncur] = tree->node[tree->node[cur[lin]].child].sibling;
			cur[lin] = tree->node[cur[lin]].child;
			ncur++;
		} else {
			for(i=lin+1; i<ncur; i++)
				cur[i-1] = cur[i];
			ncur--;
		}
		time = next;
	}
	if(fos->size>0) {
		fos->sizeBuf = fos->size;
		fos->fossilList = (TypeFossilList*) realloc((void*)fos->fossilList, fos->sizeBuf*sizeof(TypeFossilList));
	} else {
		fos->sizeBuf = 0;
		free((void*)fos->fossilList);
		fos->fossilList = NULL;
	}
	return fos;
}

/*simulate a random tree with specified birth and death rates and fossil finds on this tree with rate "fossil"*/
TypeTree *simulFossilTree(gsl_rng *rgen, double birth, double death, double fossil, double maxTime) {
	TypeTree *tree;
    int *cur, ncur, i;
	double time = 0.;
	TypeFossilFeature *fos;
	if((cur = (int*) malloc((MAX_CURRENT+1)*sizeof(int))) == NULL)
		return NULL;
	tree = (TypeTree*) malloc(sizeof(TypeTree));
	tree->sizeBuf = INC_SIZE;
	tree->node = (TypeNode*) malloc(tree->sizeBuf*sizeof(TypeNode));
	tree->time = (double*) malloc(tree->sizeBuf*sizeof(double));
	tree->maxTime = maxTime;
	tree->minTime = 0.;
	tree->size = 1;
	tree->root = 0;
	tree->time[0] = INFTY;
	tree->node[0].child = -1;
	tree->node[0].sibling = -1;
	fos = (TypeFossilFeature*) malloc(sizeof(TypeFossilFeature));
	fos->sizeBuf = INC_FOSSIL_ITEM;
	fos->fossil = (int*) malloc(tree->sizeBuf*sizeof(int));
	fos->fossil[0] = -1;
	fos->fossilList = (TypeFossilList*) malloc(fos->sizeBuf*sizeof(TypeFossilList));
	fos->size = 0;

	cur[0] = 0;
	ncur = 1;
	while(time<maxTime && ncur>0) {
		int type = getType(rgen, birth, death, fossil);
		int which = gsl_rng_uniform_int(rgen, ncur);
		double wait = gsl_ran_exponential(rgen, 1./(ncur*(birth+death+fossil)));
		time += wait;
		if(time < maxTime) {
			switch(type) {
				case 'b':
					if(ncur > MAX_CURRENT) {
						printf("too much lineages generated during simulations (%d - max %d)", ncur, MAX_CURRENT);
						freeTree(tree);
						return NULL;
//						exit(EXIT_FAILURE);
					}
					tree->time[cur[which]] = time;
					if((tree->size+1)>=tree->sizeBuf) {
						tree->sizeBuf += INC_SIZE;
						tree->node = (TypeNode*) realloc((void*)tree->node, tree->sizeBuf*sizeof(TypeNode));
						tree->time = (double*) realloc((void*)tree->time, tree->sizeBuf*sizeof(double));
						fos->fossil = (int*) realloc((void*)fos->fossil, tree->sizeBuf*sizeof(int));
					}
					tree->node[cur[which]].child = tree->size;
					tree->time[tree->size] = INFTY;
					tree->node[tree->size].child = -1;
					tree->node[tree->size].sibling = tree->size+1;
					fos->fossil[tree->size] = -1;
					tree->time[tree->size+1] = INFTY;
					tree->node[tree->size+1].child = -1;
					tree->node[tree->size+1].sibling = -1;
					fos->fossil[tree->size+1] = -1;
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
				case 'f':
					if(fos->size>=fos->sizeBuf) {
						fos->sizeBuf += INC_SIZE;
						fos->fossilList = (TypeFossilList*) realloc((void*)fos->fossilList, fos->sizeBuf*sizeof(TypeFossilList));
					}
					fos->fossilList[fos->size].time = time;
					fos->fossilList[fos->size].prec = fos->fossil[cur[which]];
					fos->fossil[cur[which]] = fos->size;
					fos->size++;
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
	if(fos->size>0) {
		fos->sizeBuf = fos->size;
		fos->fossilList = (TypeFossilList*) realloc((void*)fos->fossilList, fos->sizeBuf*sizeof(TypeFossilList));
	} else {
		fos->sizeBuf = 0;
		free((void*)fos->fossilList);
		fos->fossilList = NULL;
	}
    tree->info = (void*) fos;
	return tree;
}
