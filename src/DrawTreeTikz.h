#ifndef DrawTreeTikzF
#define DrawTreeTikzF
#include "Tree.h"
#include "DrawTreeGeneric.h"

#ifdef __cplusplus
extern "C" {
#endif

void drawTextTikz(double x0, double y0, char *text, char *mod, TypeParamDrawTreeGeneric *param);
void drawTextAngleTikz(double x0, double y0, double a, char *text, char *mod, TypeParamDrawTreeGeneric *param);
void drawDottedLineTikz(double x1, double y1, double x2, double y2, TypeParamDrawTreeGeneric *param);
void drawLineTikz(double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param);
void drawLineDotTikz(TypeRGB rgb, double alpha, double radius, double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param);
void setFunctTikz(TypeFunctDrawTreeGeneric *funct);
void startTikz(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param);
void endTikz(TypeParamDrawTreeGeneric *param);
char *sprintRGBTikz(char *buffer, TypeRGB rgb);
double getMaxLeafLabelWidthTikz(TypeTree *tree);
void startTikzStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param);
void startTikz(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param);
void endTikz(TypeParamDrawTreeGeneric *param);
void drawTextTikz(double x0, double y0, char *text, char *mod, TypeParamDrawTreeGeneric *param);
void drawLineTikz(double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param);
void fillWedgeTikz(TypeRGB rgb, double x, double y, double a, double b, TypeParamDrawTreeGeneric *param);
void drawWedgeTikz(double x, double y, double a, double b, TypeParamDrawTreeGeneric *param);
void fillGradientTikz(TypeRGB rgb0, TypeRGB rgb1, double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param);
void fillPolygonTikz(TypeRGB rgb, double alpha, double *x, double *y, int size, TypeParamDrawTreeGeneric *param);
void drawDotTikz(TypeRGB rgb, double alpha, double radius, double x, double y, TypeParamDrawTreeGeneric *param);

#ifdef __cplusplus
}
#endif

#endif
