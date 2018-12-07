#ifndef MinimizeNLOptF
#define MinimizeNLOptF

#include "Model.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct NLOPT_OPTION {
    double infSpe, supSpe, infExt, supExt, infFos, supFos, tolOptim;
    int trials, maxIter;
} TypeNLOptOption;

typedef struct ESTIMATION {
    TypeModelParam param;
    double logLikelihood;
} TypeEstimation;

int minimizeGridTwoPar(TypeLikelihoodTreeFosFunction *f, TypeTree *tree, TypeFossilFeature *fos, TypeNLOptOption *option, TypeEstimation *estim);
int minimizeBirthDeathFossilFromTreeFossil(TypeLikelihoodTreeFosFunction *f, TypeTree *tree, TypeFossilFeature *fos, TypeNLOptOption *option, TypeEstimation *estim);
int minimizeBirthDeathFromTreeFossil(TypeLikelihoodTreeFosFunction *f, TypeTree *tree, TypeFossilFeature *fos, TypeNLOptOption *option, TypeEstimation *estim);
int minimizeBirthDeathFossilFromEventList(TypeLikelihoodEventListFunction *f, TypeListEvent *event, TypeNLOptOption *option, TypeEstimation *estim);
int minimizeBirthDeathFromEventList(TypeLikelihoodEventListFunction f, TypeListEvent *event, TypeNLOptOption *option, TypeEstimation *estim);
void fprintNLoptOption(FILE *f, TypeNLOptOption *option);
void sprintNLoptOption(char *buffer, TypeNLOptOption *option);
void fprintNLoptOptionTag(FILE *f, TypeNLOptOption *option);
void fscanNLoptOptionTag(FILE *f, TypeNLOptOption *option);
#ifdef __cplusplus
}
#endif

#endif
