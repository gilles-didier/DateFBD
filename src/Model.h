#ifndef ModelF
#define ModelF

#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include "Tree.h"
#include "Fossil.h"

#define RINFTY 1E99

typedef struct MODEL_PARAM {
	double birth, death, fossil, sampling;
} TypeModelParam;

typedef struct EXTENDED_MODEL_PARAM {
    double birth, death, fossil, sampling;
    double alpha, beta, omega, delta;
    double ab, bma, bmab, amab;
    double bs, as;
} TypeExtendedModelParam;

typedef struct PIECEWISE_MODEL_PARAM {
	int size;
	double *startTime;
	TypeModelParam *param;
} TypePiecewiseModelParam;

typedef struct EVENT {
	char type;
	int n;
	double time;
} TypeEvent;

typedef struct LISTEEVENT {
	TypeEvent *list;
	int size;
	double minTime, maxTime;
} TypeListEvent;

typedef double TypeLikelihoodTreeFosFunction(TypeTree *, TypeFossilFeature *, TypeModelParam *);
typedef double TypeLikelihoodEventListFunction(TypeListEvent *, TypeModelParam *);

#ifdef __cplusplus
extern "C" {
#endif

double getLogProbExtinctClade(int *clade, double *time, double maxTime, TypeExtendedModelParam *pext);
double getLogDensClade(double val, int *clade, double *time, double maxTime, TypeExtendedModelParam *pext);
double getLogProbClade(double val, int *clade, double *time, double maxTime, TypeExtendedModelParam *pext);
double getQuantileClade(double q, double tol, int *clade, double *time, double maxTime, TypeExtendedModelParam *pext);
double getQuantile(double q, double tol, double startTime, double maxTime, TypeExtendedModelParam *param);
TypeExtendedModelParam getExtendedModelParam(TypeModelParam *param);
double getLogProbNotObservable(double t, double maxTime, TypeExtendedModelParam *param);
/*return the probability that a lineage alive at last startTime (e.g. last fossil age) goes extinct a time t without  leaving any fossil between startTime and  t*/
double getLogProbExtinct(double t, double startTime, TypeExtendedModelParam *param);
/*return the probability that a lineage alive at last startTime (e.g. last fossil age) goes extinct a time t without  leaving any fossil between startTime and  t, conditionned on the that lineage goes extinct before maxTime*/
double getLogProbExtinctCond(double t, double startTime, double maxTime, TypeExtendedModelParam *param);
double getLogDensExtinctCond(double t, double startTime, double maxTime, TypeExtendedModelParam *param);
TypePiecewiseModelParam simple2piecewise(TypeModelParam *param, double startTime, double endTime);
int getPieceIndex(double v, TypePiecewiseModelParam *param);
void printPiecewiseModel(FILE *f, TypePiecewiseModelParam *param);
void getStatEvent(TypeListEvent *event, int *b,  int *d,  int *f);
void freeListEvent(TypeListEvent *event);
void fprintEvent(FILE *f, TypeEvent *event);
void fprintListEvent(FILE *f, TypeListEvent *event);
/*get the total time*/
double getTotalTimeEvent(TypeListEvent *event);
TypeListEvent *getEventSequenceFBD(TypeTree *tree, TypeFossilFeature *fos);
TypeListEvent *getEventSequenceBD(TypeTree *tree);

double logProbEventBD(TypeListEvent *event, TypeModelParam *param);
double logProbEventFBD(TypeListEvent *event, TypeModelParam *param);


double getLogLikelihoodStadler(TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param);
double getLogLikelihoodDidier(TypeTree *tree, TypeFossilFeature *fos, TypeModelParam *param);

#ifdef __cplusplus
}
#endif

#endif
