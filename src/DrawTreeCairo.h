#ifndef DrawTreeCairoF
#define DrawTreeCairoF
#include "Tree.h"
#include "DrawTreeGeneric.h"

void drawTextCairo(double x, double y, char *text, char *al, TypeParamDrawTreeGeneric *param);
void drawTextAngleCairo(double x, double y, double a, char *text, char *al, TypeParamDrawTreeGeneric *param);
void drawDottedLineCairo(double x1, double y1, double x2, double y2, TypeParamDrawTreeGeneric *param);
void drawLineCairo(double x1, double y1, double x2, double y2, TypeParamDrawTreeGeneric *param);
void drawLineDotCairo(TypeRGB rgb, double alpha, double radius, double x1, double y1, double x2, double y2, TypeParamDrawTreeGeneric *param);
void drawDotCairo(TypeRGB rgb, double alpha, double radius, double x, double y, TypeParamDrawTreeGeneric *param);
void fillPolygonCairo(TypeRGB rgb, double alpha, double *x, double *y, int size, TypeParamDrawTreeGeneric *param);
void fillWedgeCairo(TypeRGB rgb, double x, double y, double a, double b, TypeParamDrawTreeGeneric *param);
void drawWedgeCairo(double x, double y, double a, double b, TypeParamDrawTreeGeneric *param);
void fillGradientCairo(TypeRGB rgb0, TypeRGB rgb1, double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param);
void setFunctSVG(TypeFunctDrawTreeGeneric *funct);
void startSVGStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param);
void startSVG(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param);
void setFunctPDF(TypeFunctDrawTreeGeneric *funct);
void startPDFStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param);
void startPDF(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param);
void setFunctPS(TypeFunctDrawTreeGeneric *funct);
void startPSStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param);
void startPS(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param);
void setFunctPNG(TypeFunctDrawTreeGeneric *funct);
void startPNGStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param);
void startPNG(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param);
void endPNG(TypeParamDrawTreeGeneric *param);
void endCairo(TypeParamDrawTreeGeneric *param);

#endif
