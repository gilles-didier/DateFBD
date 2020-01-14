#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include "Model.h"
#include "TreeExtras.h"

typedef struct MODEL_COEFF_BD {
    double birth, death;
    double bd;
} TypeModelCoeffBD;

typedef struct MODEL_COEFF_FBD {
    double birth, death, fossil;
    double alpha, beta, omega;
    double ab, ba, fba, bfba, afbb;
} TypeModelCoeffFBD;

static TypeModelCoeffFBD getModelCoeffFBD(TypeModelParam *param);
static double logProbObsFBD(double t, double maxTime, TypeModelCoeffFBD *coeff);
static double logProbNotObsFBD(double t, double maxTime, TypeModelCoeffFBD *coeff);
static double logDensWaitFBD(double s, double e, double maxTime, int n, TypeModelCoeffFBD *coeff);
static double logDistWaitFBD(double s, double e, double maxTime, int n, TypeModelCoeffFBD *coeff);
static double logProbBirthFBD(double t, double maxTime, TypeModelCoeffFBD *coeff);
static double logProbDeathFBD(double t, double maxTime, TypeModelCoeffFBD *coeff);
static double logProbFossilFBD(double t, double maxTime, TypeModelCoeffFBD *coeff);
static double logDistWaitOneFBD(double b, double e, double maxTime, TypeModelCoeffFBD *coeff);
static double logDensWaitOneFBD(double b, double e, double maxTime, TypeModelCoeffFBD *coeff);
static TypeModelCoeffBD getModelCoeffBD(TypeModelParam *param);
static double logProbNotObsBD(double t, double maxTime, TypeModelCoeffBD *coeff);
static double logDensWaitBD(double s, double e, double maxTime, int n, TypeModelCoeffBD *coeff);
static double logDistWaitBD(double s, double e, double maxTime, int n, TypeModelCoeffBD *coeff);
static double log_q_Stadler(double t, double c1, double c2, TypeModelParam *param);
static double log_p0_Stadler(double t, double c1, double c2, TypeModelParam *param);
//static int fillSplitTreeFossilModel(TypeTree *tcur, int n, TypeTree *tree, TypeFossilTab *ftab, TypeTree **treeList, int *size);
//static void splitTreeFossilModel(TypeTree *tree, TypeFossilFeature *fos, TypeTree ***treeList, int *size);


TypePiecewiseModelParam simple2piecewise(TypeModelParam *param, double startTime, double endTime) {
	TypePiecewiseModelParam res;
	res.size = 1;
	res.startTime = (double*) malloc(2*sizeof(double));
	res.startTime[0] = startTime;
	res.startTime[1] = endTime;
	res.param = (TypeModelParam*) malloc(sizeof(TypeModelParam));
	res.param[0] = *param;
	return res;
}

int getPieceIndex(double v, TypePiecewiseModelParam *param) {
	int a = 0, b = param->size, c;
	if(v<param->startTime[0]) {
		fprintf(stderr, "Error in 'getPieceIndex': value %.2lf too small (start time = %.2lf)\n", v, param->startTime[0]);
		return 0;
	}
	if(v>param->startTime[param->size]) {
		fprintf(stderr, "Error in 'getPieceIndex': value %.2lf too high (end time = %.2lf)\n", v, param->startTime[param->size]);
		return param->size-1;
	}
	while(b > a+1) {
		c = (a+b)/2;
		if(param->startTime[c] < v)
			a = c;
		else
			b = c;
	}
	return a;
}
			

void printPiecewiseModel(FILE *f, TypePiecewiseModelParam *param) {
	int i;
	fprintf(f, "size %d\n", param->size);
	for(i=0; i<param->size; i++)
		fprintf(f, "%lf\n%lf %lf %lf %lf\n", param->startTime[i], param->param[i].birth, param->param[i].death, param->param[i].fossil, param->param[i].sampling);
	fprintf(f, "%lf\n", param->startTime[param->size]);
}

TypeModelCoeffFBD getModelCoeffFBD(TypeModelParam *param) {
	TypeModelCoeffFBD res;
	res.birth = param->birth;
	res.death = param->death;
	res.fossil = param->fossil;
	res.omega = sqrt(pow(res.birth+res.death+res.fossil, 2.)-4.*res.birth*res.death);
	res.alpha = (res.birth+res.death+res.fossil-res.omega)/(2.*res.birth);
	res.beta = (res.birth+res.death+res.fossil+res.omega)/(2.*res.birth);
	res.ab = res.alpha*(1.-res.beta);
	res.ba = res.beta*(1.-res.alpha);
	res.fba = res.fossil+res.birth*(1.-res.alpha);
	res.bfba = res.beta*res.fba;
	res.afbb = res.alpha*(res.fossil+res.birth*(1.-res.beta));
	return res;
}

double logProbObsFBD(double t, double maxTime, TypeModelCoeffFBD *coeff) {
    return log(coeff->ba-coeff->ab*exp(-coeff->omega*(maxTime-t)))
		-log(coeff->beta-coeff->alpha*exp(-coeff->omega*(maxTime-t)));
}

double logProbNotObsFBD(double t, double maxTime, TypeModelCoeffFBD *coeff) {
    return log(coeff->alpha)+log(coeff->beta)
		+log(1-exp(-coeff->omega*(maxTime-t)))
		-log(coeff->beta-coeff->alpha*exp(-coeff->omega*(maxTime-t)));
}

double logDensWaitFBD(double s, double e, double maxTime, int n, TypeModelCoeffFBD *coeff) {
	return log((double)n)
		+((double)n-1)*log(coeff->beta-coeff->alpha*exp(-coeff->omega*(maxTime-e)))
		-((double)n)*(coeff->fba*(e-s)+log(coeff->beta-coeff->alpha*exp(-coeff->omega*(maxTime-s))))
		+log(coeff->bfba-coeff->afbb*exp(-coeff->omega*(maxTime-e)));
}

double logDistWaitFBD(double s, double e, double maxTime, int n, TypeModelCoeffFBD *coeff) {
	return ((double)n)*(-coeff->fba*(e-s)+log(coeff->beta-coeff->alpha*exp(-coeff->omega*(maxTime-e)))-log(coeff->beta-coeff->alpha*exp(-coeff->omega*(maxTime-s))));
}
double logProbBirthFBD(double t, double maxTime, TypeModelCoeffFBD *coeff) {
    return log(coeff->birth)+logProbObsFBD(t, maxTime, coeff)-log(coeff->fossil+coeff->birth*exp(logProbObsFBD(t, maxTime, coeff)));
}

double logProbDeathFBD(double t, double maxTime, TypeModelCoeffFBD *coeff) {
    return log(coeff->fossil)+logProbNotObsFBD(t, maxTime, coeff)-log(coeff->fossil+coeff->birth*exp(logProbObsFBD(t, maxTime, coeff)));
}

double logProbFossilFBD(double t, double maxTime, TypeModelCoeffFBD *coeff) {
    return log(coeff->fossil)+logProbObsFBD(t, maxTime, coeff)-log(coeff->fossil+coeff->birth*exp(logProbObsFBD(t, maxTime, coeff)));
}

double getWaitFossilProbLog(int n, double t, double x, double birth, double death, double fossil, double alpha, double beta) {
    return log(n)-((double)n)*(birth*(1.-alpha)+fossil)*x+((double)n-1.)*log(beta-alpha*exp(-birth*(beta-alpha)*(t-x)))-((double)n)*log(beta-alpha*exp(-birth*(beta-alpha)*t))+log(beta*(birth*(1.-alpha)+fossil)-alpha*(birth*(1.-beta)+fossil)*exp(-birth*(beta-alpha)*(t-x)));
}

double getEndFossilProbLog(int n, double t, double birth, double death, double fossil, double alpha, double beta) {
    return ((double)n)*(-(birth*(1.-alpha)+fossil)*t+log(beta-alpha)-log(beta-alpha*exp(-birth*(beta-alpha)*t)));
}
double logProbEventFBD(TypeListEvent *event, TypeModelParam *param) {
    double res = 0., tprec;
    int i;
    TypeModelCoeffFBD coeff;
    coeff = getModelCoeffFBD(param);
    tprec = event->minTime;
    for(i=0; i<event->size; i++) {
        res += logDensWaitFBD(tprec, event->list[i].time, event->maxTime, event->list[i].n, &coeff);
        tprec = event->list[i].time;
        switch(event->list[i].type) {
            case 'b':
                res += logProbBirthFBD(event->list[i].time, event->maxTime, &coeff);
                break;
            case 'd':
                res += logProbDeathFBD(event->list[i].time, event->maxTime, &coeff);
                break;
            case 'f':
                res += logProbFossilFBD(event->list[i].time, event->maxTime, &coeff);
                break;
        }
    }
    i = event->list[event->size-1].n;
    if(event->list[event->size-1].type == 'b')
        i++;
    if(event->list[event->size-1].type == 'd')
        i--;
    return res+logDistWaitFBD(event->list[event->size-1].time, event->maxTime, event->maxTime, i,  &coeff);
}

TypeModelCoeffBD getModelCoeffBD(TypeModelParam *param) {
	TypeModelCoeffBD res;
	res.birth = param->birth;
	res.death = param->death;
	res.bd = param->birth-param->death;
	return res;
}

double logProbNotObsBD(double t, double maxTime, TypeModelCoeffBD *coeff) {
    return log(coeff->death)
		+log(1.-exp(-coeff->bd*(maxTime-t)))
		-log(coeff->birth-coeff->death*exp(-coeff->bd*(maxTime-t)));
}

double logDensWaitBD(double s, double e, double maxTime, int n, TypeModelCoeffBD *coeff) {
	return log((double)n)
		+log(coeff->birth)
		+log(coeff->bd)
		-((double)n)*coeff->bd*(e-s)
		+((double)n-1)*log(coeff->birth-coeff->death*exp(-coeff->bd*(maxTime-e)))
		-((double)n)*log(coeff->birth-coeff->death*exp(-coeff->bd*(maxTime-s)));
}

double logDistWaitBD(double s, double e, double maxTime, int n, TypeModelCoeffBD *coeff) {
	return -((double)n)*coeff->bd*(e-s)
		+((double)n)*(log(coeff->birth-coeff->death*exp(-coeff->bd*(maxTime-e)))-log(coeff->birth-coeff->death*exp(-coeff->bd*(maxTime-s))));
}

double logProbEventBD(TypeListEvent *event, TypeModelParam *param) {
    double res = 0., tprec;
    int i;
    TypeModelCoeffBD coeff;
    coeff = getModelCoeffBD(param);
    if(event->size == 0)
        return logProbNotObsBD(event->minTime, event->maxTime, &coeff);
    tprec = event->minTime;
    for(i=0; i<event->size; i++) {
        res += logDensWaitBD(tprec, event->list[i].time, event->maxTime, event->list[i].n, &coeff);
        tprec = event->list[i].time;
    }
    return res+logDistWaitBD(event->list[event->size-1].time, event->maxTime, event->maxTime, event->list[event->size-1].n+1, &coeff);
}

/*return the sequence of speciation/extinction/fossil find events occurring in "tree", chronologically ordered in ascending order*/
/*event.list[i].n refers to the number of lineages alive just before the ith event occurs*/
TypeListEvent *getEventSequenceFBD(TypeTree *tree, TypeFossilFeature *fos) {
	double time;
	int cur[MAX_CURRENT], ncur, i, ind = 0, nbuf = INC_SIZE;
	TypeListEvent *res;
	res = (TypeListEvent*) malloc(sizeof(TypeListEvent));
	res->size = 0;
	if(tree->size == 0) {
		res->list = NULL;
		return res;
	}
	res->list = (TypeEvent*) malloc(nbuf*sizeof(TypeEvent));
	if(tree->time[tree->root] == 0.) {
		ncur = 2;
		cur[0] = tree->node[tree->root].child;
		cur[1] = tree->node[tree->node[tree->root].child].sibling;
	} else {
		ncur = 1;
		cur[0] = tree->root;
	}
	time = 0.;
	while((tree->maxTime-time)>EPSILON) {
		double tf, tbd;
		int which;
		if(ind<fos->size)
			tf = fos->fossilList[ind].time;
		else
			tf = POS_INFTY;
		tbd = tree->maxTime;
		for(i=0; i<ncur; i++) {
			if(tree->time[cur[i]] < tbd) {
				tbd = tree->time[cur[i]];
				which = i;
			}
		}
		time = utils_MIN(tf,tbd);
		if((tree->maxTime-time)>EPSILON) {
			if(res->size>=nbuf) {
				nbuf += INC_SIZE;
				res->list = (TypeEvent*) realloc((void*)res->list, nbuf*sizeof(TypeEvent));
			}
			if(tf<tbd) {
				res->list[res->size].time = tf;
				res->list[res->size].n = ncur;
				res->list[res->size].type = 'f';
				(res->size)++;
				ind++;
			} else {
				if(tf==tbd)
					ind++;
				res->list[res->size].time = tbd;
				res->list[res->size].n = ncur;
				if(tree->node[cur[which]].child>=0 && tree->node[tree->node[cur[which]].child].sibling>=0) {
					if(ncur >= MAX_CURRENT) {
						printf("too much lineages in counting events");
						exit(EXIT_FAILURE);
					}
					res->list[res->size].type = 'b';
					(res->size)++;
					cur[ncur] = tree->node[tree->node[cur[which]].child].sibling;
					cur[which] = tree->node[cur[which]].child;
					ncur++;
				} else {
					res->list[res->size].type = 'd';
					(res->size)++;
					for(i=which+1; i<ncur; i++)
						cur[i-1] = cur[i];
					ncur--;
				}
			}
		}
	}
	if(res->size)
		res->list = (TypeEvent*) realloc((void*)res->list, (res->size)*sizeof(TypeEvent));
	else {
		free((void*)res->list);
		res->list = NULL;
	}
	res->minTime = tree->minTime;
	res->maxTime = tree->maxTime;
	return res;
}

TypeListEvent *getEventSequenceBD(TypeTree *tree) {
	double time;
	int cur[MAX_CURRENT], ncur, i, nbuf = INC_SIZE;
	TypeListEvent *res;

	res = (TypeListEvent*) malloc(sizeof(TypeListEvent));
	res->size = 0;
	if(tree == NULL || tree->size == 0) {
		fprintf(stderr, "Empty tree\n");
		exit(1);
		res->list = NULL;
		return res;
	}
	res->list = (TypeEvent*) malloc(nbuf*sizeof(TypeEvent));
	if(tree->time[tree->root] == tree->minTime) {
		if(tree->node[tree->root].child != NOSUCH) {
			if(res->size>=nbuf) {
				nbuf += INC_SIZE;
				res->list = (TypeEvent*) realloc((void*)res->list, nbuf*sizeof(TypeEvent));
			}
			res->list[res->size].time = tree->minTime;
			res->list[res->size].n = 1;
			res->list[res->size].type = 'b';
			(res->size)++;
			ncur = 2;
			cur[0] = tree->node[tree->root].child;
			cur[1] = tree->node[tree->node[tree->root].child].sibling;
		} else {
			if(res->size>=nbuf) {
				nbuf += INC_SIZE;
				res->list = (TypeEvent*) realloc((void*)res->list, nbuf*sizeof(TypeEvent));
			}
			res->list[res->size].time = tree->minTime;
			res->list[res->size].n = 1;
			res->list[res->size].type = 'd';
			(res->size)++;
		}
	} else {
		ncur = 1;
		cur[0] = tree->root;
	}
	time = 0;
	while((tree->maxTime-time)>EPSILON) {
		int which;
		time = tree->maxTime;
		for(i=0; i<ncur; i++) {
			if(tree->time[cur[i]] < time) {
				time = tree->time[cur[i]];
				which = i;
			}
		}
		if((tree->maxTime-time)>EPSILON) {
			if(res->size>=nbuf) {
				nbuf += INC_SIZE;
				res->list = (TypeEvent*) realloc((void*)res->list, nbuf*sizeof(TypeEvent));
			}
			res->list[res->size].time = time;
			res->list[res->size].n = ncur;
			if(tree->node[cur[which]].child != NOSUCH) {
				if(tree->node[tree->node[cur[which]].child].sibling != NOSUCH) {
					if(ncur >= MAX_CURRENT) {
						printf("too much lineages in counting events");
						exit(EXIT_FAILURE);
					}
					res->list[res->size].type = 'b';
					(res->size)++;
					cur[ncur] = tree->node[tree->node[cur[which]].child].sibling;
					cur[which] = tree->node[cur[which]].child;
					ncur++;
				} else {
					fprintf(stderr, "Execution error: node %d with a single child\n", cur[which]);
					exit(1);
				}
			} else {
				res->list[res->size].type = 'd';
				(res->size)++;
				for(i=which+1; i<ncur; i++)
					cur[i-1] = cur[i];
				ncur--;
			}
		}
	}
	if(res->size)
		res->list = (TypeEvent*) realloc((void*)res->list, (res->size)*sizeof(TypeEvent));
	else {
		free((void*)res->list);
		res->list = NULL;
	}
	res->minTime = tree->minTime;
	res->maxTime = tree->maxTime;
	return res;
}

/*get the total time*/
double getTotalTimeEvent(TypeListEvent *event) {
    double tot = 0., tprec = 0.;
    int n;
    for(n=0; n<event->size; n++) {
        tot += (event->list[n].time-tprec)*event->list[n].n;
        tprec = event->list[n].time;
    }
    if(event->size>0) {
        n = event->list[event->size-1].n;
        if(event->list[event->size-1].type == 'b')
            n++;
        if(event->list[event->size-1].type == 'd')
            n--;
    } else
        n = 1;
    return tot+n*(event->maxTime-tprec);
}


/*shift all the times to set first event times to 0*/
void offsetTimeEvents(TypeListEvent *event) {
    int i;
    double zero;

    if(event->size ==0)
        return;
    zero = event->list[0].time;
    for(i=0; i<event->size; i++)
        event->list[i].time -= zero;
    event->maxTime -= zero;
}

/*set the number of events variables b, d, f for each type*/
void getStatEvent(TypeListEvent *event, int *b,  int *d,  int *f) {
    int n;

    *b = 0;
    *d = 0;
    *f = 0;
    for(n=0; n<event->size; n++) {
        switch(event->list[n].type) {
            case 'b':
                (*b)++;
                break;
            case 'd':
                (*d)++;
                break;
            case 'f':
                (*f)++;
                break;
            default:
            ;
        }
    }
}

/*free list of events*/
void freeListEvent(TypeListEvent *event) {
    if(event->size)
        free((void*)event->list);
    free((void*)event);
}

/*print event*/
void fprintEvent(FILE *f, TypeEvent *event) {
    fprintf(f, "%c\t%d\t%.2lE", event->type, event->n, event->time);
}

/*print list of events*/
void fprintListEvent(FILE *f, TypeListEvent *event) {
    int i;
    for(i=0; i<event->size; i++) {
        fprintEvent(f, &(event->list[i]));
        fprintf(f, "\n");
    }
}

double log_q_Stadler(double t, double c1, double c2, TypeModelParam *param) {
	return log(2.*(1-pow(c2, 2.))+exp(-c1*t)*pow(1.-c2, 2.)+exp(c1*t)*pow(1.+c2, 2.));
}

double log_p0_Stadler(double t, double c1, double c2, TypeModelParam *param) {
	return log(param->birth+param->death+param->fossil+(c1*(exp(-c1*t)*(1-c2)-(1+c2)))/(exp(-c1*t)*(1-c2)+(1+c2)))-log(2.)-log(param->birth);
}

double log_p1_Stadler(double t, double c1, double c2, TypeModelParam *param) {
	return log(4.)-log_q_Stadler(t, c1, c2, param);
}

double getLogLikelihoodStadler(TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param) {
	double *x, *y, c1, c2, ll;
	int i, n, m, e;
	x = (double*) malloc(((tree->size+1)/2)*sizeof(double));
	y = (double*) malloc(((tree->size+1)/2)*sizeof(double));
	n=0; m=0; e=0;
	for(i=0; i<tree->size; i++) {
		if(tree->node[i].child != NOSUCH) {
			x[n++] = tree->maxTime-tree->time[i];
		} else {
			if(fabs(tree->maxTime-tree->time[i])>EPSILON) {
				y[m++] = tree->maxTime-tree->time[i];
			} else
				e++;
		}
	}
	c1 = sqrt(pow(param->birth-param->death-param->fossil, 2.)+4.*param->birth*param->fossil),
	c2 = (param->birth+param->death+param->fossil)/c1;
	ll = ((double)n)*log(param->birth)
		+((double)fos->size)*log(param->fossil)
		+((double)e)*log(4.)
		-log_q_Stadler(tree->minTime, c1, c2, param);
	for(i=0; i<n; i++)
		ll -= log_q_Stadler(x[i], c1, c2, param);
	for(i=0; i<m; i++)
		ll += log_p0_Stadler(y[i], c1, c2, param)+log_q_Stadler(y[i], c1, c2, param);
	free((void*)x);
	free((void*)y);
	return ll;
}

double logDistWaitOneFBD(double b, double e, double maxTime, TypeModelCoeffFBD *coeff) {
	return -coeff->fba*(e-b)
		+log(coeff->beta-coeff->alpha*exp(-coeff->omega*(maxTime-e)))
		-log(coeff->beta-coeff->alpha*exp(-coeff->omega*(maxTime-b)));
}

/*return the log of the wait density divided by coeff->beta*(coeff->fossil+coeff->birth*(1-coeff->alpha))-coeff->alpha*(coeff->fossil+coeff->birth*(1-coeff->beta))*exp(-coeff->omega*(maxTime-e))*/
double logDensWaitOneFBD(double b, double e, double maxTime, TypeModelCoeffFBD *coeff) {
	return -coeff->fba*(e-b)
		-log(coeff->beta-coeff->alpha*exp(-coeff->omega*(maxTime-b)));
}

/* x/n times/number of divergences
 * y/m times/number of internal fossils
 * z/k times/number of external fossils
 */
double getLogLikelihoodDidier(TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param) {
	double *x, *y, *z, ll;
	int i, j, n, m, k;
	TypeModelCoeffFBD coeff;
    //TypeTree **treeList;
    //int size;
    	
	coeff = getModelCoeffFBD(param);
	x = (double*) malloc(((tree->size+1)/2)*sizeof(double));
	y = (double*) malloc(((tree->size+1)/2)*sizeof(double));
	z = (double*) malloc(fos->size*sizeof(double));
	n=0; m=0; k=0;
	for(i=0; i<tree->size; i++) {
		if(tree->node[i].child != NOSUCH) {
			x[n++] = tree->time[i];
		} else {
			if(fabs(tree->maxTime-tree->time[i])>EPSILON) {
				y[m++] = tree->time[i];
			}
		}
		for(j=fos->fossil[i]; j!=NOSUCH; j = fos->fossilList[j].prec)
			if(fabs(fos->fossilList[j].time-tree->time[i])>EPSILON)
				z[k++] = fos->fossilList[j].time;
	}
	if(tree->parent == NULL)
		tree->parent = getParent(tree);
	ll = 0.;
	for(i=0; i<tree->size; i++) {
		double start = (tree->parent[i] != NOSUCH)?tree->time[tree->parent[i]]:tree->minTime, end = tree->time[i];
		j = fos->fossil[i];
		if(j!=NOSUCH && fabs(fos->fossilList[j].time-tree->time[i])<=EPSILON)
			j = fos->fossilList[j].prec;
		for(; j!=NOSUCH; j = fos->fossilList[j].prec) {
			if(fabs(end-tree->maxTime)<=EPSILON)
				ll += logDistWaitOneFBD(fos->fossilList[j].time, end, tree->maxTime, &coeff);
			else
				ll += logDensWaitOneFBD(fos->fossilList[j].time, end, tree->maxTime, &coeff);
			end = fos->fossilList[j].time;
		}
		if(fabs(end-tree->maxTime)<=EPSILON)
			ll += logDistWaitOneFBD(start, end, tree->maxTime, &coeff);
		else
			ll += logDensWaitOneFBD(start, end, tree->maxTime, &coeff);
	}
	for(i=0; i<n; i++)
		ll += log(coeff.birth)+log(coeff.ba-coeff.ab*exp(-coeff.omega*(tree->maxTime-x[i])));
	for(i=0; i<m; i++)
		ll += log(coeff.fossil)+log(coeff.alpha)+log(coeff.beta)+log(1.-exp(-coeff.omega*(tree->maxTime-y[i])));
	for(i=0; i<k; i++)
		ll += log(coeff.fossil)+log(coeff.ba-coeff.ab*exp(-coeff.omega*(tree->maxTime-z[i])));
	free((void*)x);
	free((void*)y);
	free((void*)z);
    //splitTreeFossilModel(tree, fos, &treeList, &size);
    //for(i=0; i<size; i++)
		//if(treeList[i] != NULL) {
			//ll -= gsl_sf_lnfact((double)(treeList[i]->size+1)/2);
			//freeTree(treeList[i]);
		//}
	//free((void*)treeList);
	return ll;
}

//int fillSplitTreeFossilModel(TypeTree *tcur, int n, TypeTree *tree, TypeFossilTab *ftab, TypeTree **treeList, int *size) {
    //int curInd = tcur->size++;
    //if(ftab[n].size>0) {
        //tcur->time[curInd] = ftab[n].time[0];
        //tcur->node[curInd].child = -1;
        //if(ftab[n].size == 1 && tree->node[n].child == NOSUCH && (tree->time[n] == ftab[n].time[0] || tree->time[n]==NO_TIME)) { /*empty tree*/
            //treeList[(*size)] = newTree(0);
            //treeList[(*size)]->root = 0;
            //treeList[(*size)]->minTime = ftab[n].time[0];
            //treeList[(*size)]->minTimeInt.inf = ftab[n].time[0];
            //treeList[(*size)]->minTimeInt.sup = ftab[n].time[0];
            //treeList[(*size)]->maxTime = tree->maxTime;
            //(*size)++;
            //ftab[n].size--;
            //ftab[n].time++;
        //} else  { /*tree not empty*/
            //treeList[(*size)] = newTree(tree->size);
            //treeList[(*size)]->root = 0;
            //treeList[(*size)]->minTime = tcur->time[curInd];
            //treeList[(*size)]->minTimeInt.inf = tcur->time[curInd];
            //treeList[(*size)]->minTimeInt.sup = tcur->time[curInd];
            //treeList[(*size)]->maxTime = tree->maxTime;
            //(*size)++;
            //ftab[n].size--;
            //ftab[n].time++;
            //fillSplitTreeFossilModel(treeList[(*size)-1], n, tree, ftab, treeList, size);
        //}
    //} else {
        //tcur->time[curInd] = tree->time[n];
        //if(tree->node[n].child!=NOSUCH) {
            //int c, prec;
            //tcur->node[curInd].child = fillSplitTreeFossilModel(tcur, tree->node[n].child, tree, ftab, treeList, size);
            //prec = tcur->node[curInd].child;
            //for(c=tree->node[tree->node[n].child].sibling; c!=NOSUCH; c=tree->node[c].sibling) {
                //tcur->node[prec].sibling = fillSplitTreeFossilModel(tcur, c, tree, ftab, treeList, size);
                //prec = tcur->node[prec].sibling;
            //}
            //tcur->node[prec].sibling = NOSUCH;
        //} else {
            //tcur->node[curInd].child = NOSUCH;
        //}
    //}
    //return curInd;
//}

///*split the tree "tree" each time a fossill occurs*/
//void splitTreeFossilModel(TypeTree *tree, TypeFossilFeature *fos, TypeTree ***treeList, int *size) {
    //int i;
    //TypeFossilTab *ftab, *ftmp;
    //ftab = listToFossilTab(fos, tree->size);
    //ftmp = (TypeFossilTab*) malloc(tree->size*sizeof(TypeFossilTab));
    //for(i=0; i<tree->size; i++) {
        //ftmp[i].size = ftab[i].size;
        //ftmp[i].time = ftab[i].time;
    //}
    //*treeList = (TypeTree **) malloc((fos->size+1)*sizeof(TypeTree*));
    //(*treeList)[0] = newTree(tree->size);
    //(*treeList)[0]->root = 0;
    //(*treeList)[0]->node[0].sibling = NOSUCH;
    //(*treeList)[0]->minTime = tree->minTime;
    //(*treeList)[0]->minTimeInt.inf = tree->minTimeInt.inf;
    //(*treeList)[0]->minTimeInt.sup = tree->minTimeInt.sup;
    //(*treeList)[0]->maxTime = tree->maxTime;
    //*size = 1;
    //fillSplitTreeFossilModel((*treeList)[0], tree->root, tree, ftmp, *treeList, size);
    //for(i=0; i<tree->size; i++)
        //if(ftab[i].time != NULL)
            //free((void*)ftab[i].time);
    //free((void*)ftab);
    //free((void*)ftmp);
    //for(i=0; i<*size; i++)
        //reallocTree((*treeList)[i]->size, (*treeList)[i]);
//}

/*
double log_q_Stadler(double t, TypeModelParam *param) {
	double c1, c2;
	c1 = sqrt(pow(param->birth-param->death-param->fossil, 2.)-4.*param->birth*param->fossil),
	c2 = -(param->birth-param->death-2.*param->sampling*param->birth-param->fossil)/c1;
	return log(2.)+ log(1-pow(c2, 2.)+exp(-c1*t)*(1+pow(c2, 2.)));
}

double log_p0_Stadler(double t, TypeModelParam *param) {
	double c1, c2;
	c1 = sqrt(pow(param->birth-param->death-param->fossil, 2.)+4.*param->birth*param->fossil),
	c2 = -(param->birth-param->death-2.*param->sampling*param->birth-param->fossil)/c1;
	return log(param->birth+param->death+param->fossil+(c1*(exp(-c1*t)*(1-c2)-(1+c2)))/(exp(-c1*t)*(1-c2)+(1+c2)))-log(2.)-log(param->birth);
}

double getLogLikelihoodStadler(TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param) {
	double *x, *y, ll;
	int i, n, m;
	x = (double*) malloc(((tree->size+1)/2)*sizeof(double));
	y = (double*) malloc(((tree->size+1)/2)*sizeof(double));
	n=0; m=0;
	for(i=0; i<tree->size; i++) {
		if(tree->node[i].child == NOSUCH) {
			if(fabs(tree->maxTime-tree->time[i])<=EPSILON) {
				x[n++] = tree->time[i];
			} else {
				y[m++] = tree->time[i];
			}
		}
	}
	ll = ((double)(n+m-1))*log(param->birth)
		+((double)fos->size)*log(param->fossil)
		+((double)n)*(log(4.)+log(param->sampling))
		-log_q_Stadler(tree->minTime, param);
	for(i=0; i<(n+m-1); i++)
		ll -= log_q_Stadler(x[i], param);
	for(i=0; i<m; i++)
		ll += log_p0_Stadler(y[i], param)+log_q_Stadler(y[i], param);
	free((void*)x);
	free((void*)y);
	return ll;
}
*/
