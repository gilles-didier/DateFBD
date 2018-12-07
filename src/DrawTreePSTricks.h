#ifndef DrawTreePstricksF
#define DrawTreePstricksF
#include "Tree.h"
#include "DrawTreeGeneric.h"

#ifdef __cplusplus
extern "C" {
#endif



void drawTextPSTricks(double x0, double y0, char *text, char *mod, TypeParamDrawTreeGeneric *param);
void drawTextAnglePSTricks(double x0, double y0, double a, char *text, char *mod, TypeParamDrawTreeGeneric *param);
void drawDottedLineCairo(double x1, double y1, double x2, double y2, TypeParamDrawTreeGeneric *param);
void drawLinePSTricks(double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param);
void setFunctPSTricks(TypeFunctDrawTreeGeneric *funct);
void startPSTricks(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param);
void endPSTricks(TypeParamDrawTreeGeneric *param);
char *sprintRGBPSTricks(char *buffer, TypeRGB rgb);
double getMaxLeafLabelWidthPSTricks(TypeTree *tree);
void drawLineDotPSTricks(TypeRGB rgb, double alpha, double radius, double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param);
void drawDotPSTricks(TypeRGB rgb, double alpha, double radius, double x, double y, TypeParamDrawTreeGeneric *param);

void startPSTricksStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param);
void startPSTricks(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param);
void endPSTricks(TypeParamDrawTreeGeneric *param);
void drawTextPSTricks(double x0, double y0, char *text, char *mod, TypeParamDrawTreeGeneric *param);
void drawLinePSTricks(double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param);
void fillWedgePSTricks(TypeRGB rgb, double x, double y, double a, double b, TypeParamDrawTreeGeneric *param);
void drawWedgePSTricks(double x, double y, double a, double b, TypeParamDrawTreeGeneric *param);
void fillGradientPSTricks(TypeRGB rgb0, TypeRGB rgb1, double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param);
void fillPolygonPSTricks(TypeRGB rgb, double alpha, double *x, double *y, int size, TypeParamDrawTreeGeneric *param);

#ifdef __cplusplus
}
#endif

#endif
