#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <nlopt.h>

#include "Utils.h"

#include "MinimizeNLOpt.h"


//#define NLOPT_ALGO NLOPT_GN_ISRES
//#define NLOPT_ALGO NLOPT_GN_ESCH
#define NLOPT_ALGO NLOPT_LN_BOBYQA
//#define NLOPT_ALGO NLOPT_LN_COBYLA
//#define NLOPT_ALGO NLOPT_AUGLAG


#define MINVAL 0.01
#define INFTY 1E99
#define RINFTY 1E99
#define DEF 10
#define MIN_VAL 0.000001
#define TOLERANCE_CONSTRAINT 0.000000001
#define TOLERANCE_OPTIM 0.001



static double toMinimizeTreeFossil(unsigned n, const double *x, double *grad, void *data);
static double toMinimizeTreeFossilFixed(unsigned n, const double *x, double *grad, void *data);
static double toMinimizeEventList(unsigned n, const double *x, double *grad, void *data);
static double toMinimizeEventListFossilFixed(unsigned n, const double *x, double *grad, void *data);
static double ratesConstraint(unsigned n, const double *x, double *grad, void *data);

typedef struct MINIMIZATION_TREE_FOS_DATA {
    TypeTree *tree;
    TypeFossilFeature *fos;
    TypeLikelihoodTreeFosFunction *likelihood;
} TypeMinimizationTreeFossilData;

typedef struct MINIMIZATION_TREE_FOS_FIXED_DATA {
    TypeTree *tree;
    TypeFossilFeature *fos;
    double fossil;
    TypeLikelihoodTreeFosFunction *likelihood;
} TypeMinimizationTreeFossilFixedData;

typedef struct MINIMIZATION_EVENT_LIST_FOS_FIXED_DATA {
    TypeListEvent *event;
    double fossil;
    TypeLikelihoodEventListFunction *likelihood;
} TypeMinimizationEventListFossilFixedData;

typedef struct MINIMIZATION_EVENT_LIST_DATA {
    TypeListEvent *event;
    TypeLikelihoodEventListFunction *likelihood;
} TypeMinimizationEventListData;

#define TAG_SPE "SPE"
#define TAG_EXT "EXT"
#define TAG_FOS "FOS"
#define TAG_TRI "TRI"
#define TAG_TOL "TOL"
#define TAG_ITE "ITE"
#define SIZE_TAG 20
#define SIZE_VAL 100

void fprintNLoptOptionTag(FILE *f, TypeNLOptOption *option) {
    fprintf(f, ":%s [%lE;%lE]\n", TAG_SPE, option->infSpe, option->supSpe);
    fprintf(f, ":%s [%lE;%lE]\n", TAG_EXT, option->infExt, option->supExt);
    fprintf(f, ":%s [%lE;%lE]\n", TAG_FOS, option->infFos, option->supFos);
    fprintf(f, ":%s %d\n", TAG_TRI, option->trials);
    fprintf(f, ":%s %lE\n", TAG_TOL, option->tolOptim);
    fprintf(f, ":%s %d\n", TAG_ITE, option->maxIter);
}

void fscanNLoptOptionTag(FILE *f, TypeNLOptOption *option) {
    char c, tag[SIZE_TAG+1], val[SIZE_VAL+1];
    for(c=fgetc(f); c!=EOF && isspace(c); c=fgetc(f));
    while(c == ':') {
        int i;
        c=fgetc(f);
        for(i=0; c!=EOF && !isspace(c) && i<SIZE_TAG; c=fgetc(f))
            tag[i++] = c;
        tag[i] = '\0';
        if(i>=SIZE_TAG) {
            fprintf(stderr, "Error when reading an optimizer options file - Tag too long:\n%s...\n", tag);
            exit(1);
        }
        for(; c!=EOF && isspace(c); c=fgetc(f));
        for(i=0; c!=EOF && !isspace(c) && i<SIZE_VAL; c=fgetc(f))
            val[i++] = c;
        val[i] = '\0';
        if(i>=SIZE_VAL) {
            fprintf(stderr, "Error when reading an optimizer options file - value too long:\n%s...\n", val);
            exit(1);
        }
        if(strcmp(tag, TAG_SPE) == 0)
            toInterval(val, &(option->infSpe), &(option->supSpe));
        if(strcmp(tag, TAG_EXT) == 0)
            toInterval(val, &(option->infExt), &(option->supExt));
        if(strcmp(tag, TAG_FOS) == 0)
            toInterval(val, &(option->infFos), &(option->supFos));
        if(strcmp(tag, TAG_TRI) == 0)
            option->trials = atoi(val);
        if(strcmp(tag, TAG_TOL) == 0)
            option->tolOptim = atof(val);
        if(strcmp(tag, TAG_ITE) == 0)
            option->maxIter = atoi(val);
        for(; c!=EOF && isspace(c); c=fgetc(f));
    }
}

void fprintNLoptOption(FILE *f, TypeNLOptOption *option) {
    fprintf(f, "Speciation rates are sampled in [%.2lE:%.2lE]\n", option->infSpe, option->supSpe);
    fprintf(f, "Extinction rates are sampled in [%.2lE:%.2lE]\n", option->infExt, option->supExt);
    fprintf(f, "Fossil rates are sampled in [%.2lE:%.2lE]\n", option->infFos, option->supFos);
    fprintf(f, "Optimizer runs %d trials and stops with tolerance %.lE or after more than %d iterations.\n", option->trials, option->tolOptim, option->maxIter);
}

void sprintNLoptOption(char *buffer, TypeNLOptOption *option) {
    buffer += sprintf(buffer, "Speciation rates are sampled in [%.2lE:%.2lE]\n", option->infSpe, option->supSpe);
    buffer += sprintf(buffer, "Extinction rates are sampled in [%.2lE:%.2lE]\n", option->infExt, option->supExt);
    buffer += sprintf(buffer, "Fossil rates are sampled in [%.2lE:%.2lE]\n", option->infFos, option->supFos);
    buffer += sprintf(buffer, "Optimizer runs %d trials and stops with tolerance %.lE or after more than %d iterations.\n", option->trials, option->tolOptim, option->maxIter);
}

double toMinimizeTreeFossil(unsigned n, const double *x, double *grad, void *data) {
    TypeModelParam param;
    param.birth = x[0];
    param.death = x[1];
    param.fossil = x[2];
    param.sampling = 1.;
    return -((TypeMinimizationTreeFossilData*)data)->likelihood(((TypeMinimizationTreeFossilData*)data)->tree, ((TypeMinimizationTreeFossilData*)data)->fos, &param);
}

double toMinimizeTreeFossilFixed(unsigned n, const double *x, double *grad, void *data) {
    TypeModelParam param;
    param.birth = x[0];
    param.death = x[1];
    param.fossil = ((TypeMinimizationTreeFossilFixedData*)data)->fossil;
    param.sampling = 1.;
    return -((TypeMinimizationTreeFossilFixedData*)data)->likelihood(((TypeMinimizationTreeFossilFixedData*)data)->tree, ((TypeMinimizationTreeFossilFixedData*)data)->fos, &param);
}

double toMinimizeEventList(unsigned n, const double *x, double *grad, void *data) {
    TypeModelParam param;
    param.birth = x[0];
    param.death = x[1];
    param.fossil = x[2];
    return -((TypeMinimizationEventListData*)data)->likelihood(((TypeMinimizationEventListData*)data)->event, &param);
}

double toMinimizeEventListFossilFixed(unsigned n, const double *x, double *grad, void *data) {
    TypeModelParam param;
    param.birth = x[0];
    param.death = x[1];
    param.fossil = ((TypeMinimizationEventListFossilFixedData*)data)->fossil;
    return -((TypeMinimizationEventListFossilFixedData*)data)->likelihood(((TypeMinimizationEventListFossilFixedData*)data)->event, &param);
}

double ratesConstraint(unsigned n, const double *x, double *grad, void *data) {
//	fprintf(stderr, "constraint %.2le\n", x[1]-x[0]);
	return x[1]-x[0];
}


int minimizeBirthDeathFossilFromTreeFossil(TypeLikelihoodTreeFosFunction *f, TypeTree *tree, TypeFossilFeature *fos, TypeNLOptOption *option, TypeEstimation *estim) {
    double lb[3], x[3], minLikelihood;
    nlopt_opt opt;
    TypeMinimizationTreeFossilData data;
    int result, t;
    data.tree = tree;
    data.fos = fos;
    data.likelihood = f;
    lb[0] = 0.;
    lb[1] = 0.;
    lb[2] = 0.;
    opt = nlopt_create(NLOPT_ALGO, 3); /* algorithm and dimensionality */
	nlopt_add_inequality_constraint(opt, ratesConstraint, NULL, 1e-8);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_min_objective(opt, toMinimizeTreeFossil, &data);
    nlopt_set_xtol_abs1(opt, option->tolOptim);
    nlopt_set_maxeval(opt, option->maxIter);
    estim->logLikelihood = INFTY;
    for(t=0; t<option->trials; t++) {
        x[1] = option->infExt+UNIF_RAND*(option->supExt-option->infExt);
        x[0] = x[1]+UNIF_RAND*(option->supSpe-x[1]);
        x[2] = option->infFos+UNIF_RAND*(option->supFos-option->infFos);
        if(((result = nlopt_optimize(opt, x, &minLikelihood)) >= 0) && minLikelihood < estim->logLikelihood) {
            estim->logLikelihood = minLikelihood;
            estim->param.birth = x[0];
            estim->param.death = x[1];
            estim->param.fossil = x[2];
        }
    }
    estim->logLikelihood = -estim->logLikelihood;
    nlopt_destroy(opt);
    return result;
}

int minimizeBirthDeathFromTreeFossil(TypeLikelihoodTreeFosFunction *f, TypeTree *tree, TypeFossilFeature *fos, TypeNLOptOption *option, TypeEstimation *estim) {
    double lb[2], x[2], minLikelihood;
    nlopt_opt opt;
    TypeMinimizationTreeFossilFixedData data;
    int result, t;
    data.tree = tree;
    data.fos = fos;
    data.fossil = estim->param.fossil;
    data.likelihood = f;
    lb[0] = 0.;
    lb[1] = 0.;
    opt = nlopt_create(NLOPT_ALGO, 2); /* algorithm and dimensionality */
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_min_objective(opt, toMinimizeTreeFossilFixed, &data);
    nlopt_set_xtol_abs1(opt, option->tolOptim);
    nlopt_set_maxeval(opt, option->maxIter);
    estim->logLikelihood = INFTY;
    for(t=0; t<option->trials; t++) {
        x[0] = option->infSpe+UNIF_RAND*(option->supSpe-option->infSpe);
        x[1] = option->infExt+UNIF_RAND*(option->supExt-option->infExt);
        if(((result = nlopt_optimize(opt, x, &minLikelihood)) >= 0) && minLikelihood < estim->logLikelihood) {
            estim->logLikelihood = minLikelihood;
            estim->param.birth = x[0];
            estim->param.death = x[1];
        }
    }
    estim->logLikelihood = -estim->logLikelihood;
    nlopt_destroy(opt);
    return result;
}

int minimizeBirthDeathFossilFromEventList(TypeLikelihoodEventListFunction *f, TypeListEvent *event, TypeNLOptOption *option, TypeEstimation *estim) {
    double x[3], lb[3], minLikelihood;
    nlopt_opt opt;
    TypeMinimizationEventListData data;
    int result, t;

    data.event = event;
    data.likelihood = f;
    lb[0] = 0.;
    lb[1] = 0.;
    lb[2] = 0.;
    opt = nlopt_create(NLOPT_ALGO, 3); /* algorithm and dimensionality */
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_min_objective(opt, toMinimizeEventList, &data);
    nlopt_set_xtol_abs1(opt, option->tolOptim);
    nlopt_set_maxeval(opt, option->maxIter);
    estim->logLikelihood = INFTY;
    for(t=0; t<option->trials; t++) {
        x[0] = option->infSpe+UNIF_RAND*(option->supSpe-option->infSpe);
        x[1] = option->infExt+UNIF_RAND*(option->supExt-option->infExt);
        x[2] = option->infFos+UNIF_RAND*(option->supFos-option->infFos);
        if(((result = nlopt_optimize(opt, x, &minLikelihood)) >= 0) && minLikelihood < estim->logLikelihood) {
            estim->logLikelihood = minLikelihood;
            estim->param.birth = x[0];
            estim->param.death = x[1];
            estim->param.fossil = x[2];
        }
    }
    estim->logLikelihood = -estim->logLikelihood;
    nlopt_destroy(opt);
    return result;
}


int minimizeBirthDeathFromEventList(TypeLikelihoodEventListFunction *f, TypeListEvent *event, TypeNLOptOption *option, TypeEstimation *estim) {
    double lb[2], x[2], minLikelihood;
    nlopt_opt opt;
    TypeMinimizationEventListFossilFixedData data;
    int result, t;

    data.event = event;
    data.likelihood = f;
    data.fossil = estim->param.fossil;
    lb[0] = 0.;
    lb[1] = 0.;
    opt = nlopt_create(NLOPT_ALGO, 2); /* algorithm and dimensionality */
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_min_objective(opt, toMinimizeEventListFossilFixed, &data);
    nlopt_set_lower_bounds1(opt, 0.);
    nlopt_set_xtol_abs1(opt, option->tolOptim);
    nlopt_set_maxeval(opt, option->maxIter);
    estim->logLikelihood = INFTY;
    for(t=0; t<option->trials; t++) {
        x[0] = option->infSpe+UNIF_RAND*(option->supSpe-option->infSpe);
        x[1] = option->infExt+UNIF_RAND*(option->supExt-option->infExt);
        if(((result = nlopt_optimize(opt, x, &minLikelihood)) >= 0) && minLikelihood < estim->logLikelihood) {
            estim->logLikelihood = minLikelihood;
            estim->param.birth = x[0];
            estim->param.death = x[1];
       }
    }
    estim->logLikelihood = -estim->logLikelihood;
    nlopt_destroy(opt);
    return result;
}
