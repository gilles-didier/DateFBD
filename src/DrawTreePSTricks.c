#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <locale.h>
#include "Utils.h"
#include "DrawTreePSTricks.h"

#define NB_MIN 6
#define OFFSET 10
#define TICK_LENGTH 4
#define FONT_NAME "Helvetica"
//#define LABEL_SEP 10.
#define MAX_STRING_SIZE 200
//#define STANDARD_WIDTH 500
#define CHAR_WIDTH 3
#define LABEL_SEP 10.
#define FONT_SIZE 9.
#define STANDARD_WIDTH 500
#define MY_PI acos(-1)



char *sprintRGBPSTricks(char *buffer, TypeRGB rgb) {
	sprintf(buffer, "{[rgb]{%.3lf %.3lf %.3lf}}", rgb.red, rgb.green, rgb.blue);
	return buffer;
}

void drawTextPSTricks(double x0, double y0, char *text, char *mod, TypeParamDrawTreeGeneric *param) {
    int i;
    char *tmp = strdpl(text);
    for(i=0; tmp[i]!='\0'; i++)
        if(tmp[i] == '_')
            tmp[i] = ' ';
    fprintf((FILE*)param->info, "\\rput[%s](%.2lf, %.2lf){\\small \\textcolor[rgb]{%.3lf %.3lf %.3lf}{%s}}\n", mod, x0, param->height-y0, param->curgb.red, param->curgb.green, param->curgb.blue, tmp);
    free((void*)tmp);
}

void drawTextAnglePSTricks(double x0, double y0, double a, char *text, char *mod, TypeParamDrawTreeGeneric *param) {
	int i;
	char *tmp = strdpl(text);
	for(i=0; tmp[i]!='\0'; i++)
		if(tmp[i] == '_')
			tmp[i] = ' ';
	fprintf((FILE*)param->info, "\\rput[%s]{%.2lf}(%.2lf,%.2lf){\\small \\textcolor[rgb]{%.3lf %.3lf %.3lf}{%s}}\n", mod, -(180.*a)/MY_PI, x0, param->height-y0, param->curgb.red, param->curgb.green, param->curgb.blue, tmp);
	free((void*)tmp);
}

void drawLinePSTricks(double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param) {
	char buffer[MAX_STRING_SIZE];
    fprintf((FILE*)param->info, "\\psline[linewidth=1pt,linecolor=%s](%.2lf, %.2lf)(%.2lf, %.2lf)\n", sprintRGBPSTricks(buffer, param->curgb), x0, param->height-y0, x1, param->height-y1);
}

void drawDottedLinePSTricks(double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param) {
	char buffer[MAX_STRING_SIZE];
    fprintf((FILE*)param->info, "\\psline[linestyle=dotted,linewidth=1pt,linecolor=%s](%.2lf, %.2lf)(%.2lf, %.2lf)\n", sprintRGBPSTricks(buffer, param->curgb), x0, param->height-y0, x1, param->height-y1);
}

void fillWedgePSTricks(TypeRGB rgb, double x, double y, double a, double b, TypeParamDrawTreeGeneric *param) {
	char buffer[MAX_STRING_SIZE];
	fprintf((FILE*)param->info, "\\pswedge[linestyle=none,fillstyle=solid,fillcolor=%s,fillcolor=%s](%lf,%lf){%lfpt}{%lf}{%lf}\n", sprintRGBPSTricks(buffer, rgb), sprintRGBPSTricks(buffer, param->curgb), x, param->height-y, param->radius, 90*a/asin(1), 90*b/asin(1));
}

void drawWedgePSTricks(double x, double y, double a, double b, TypeParamDrawTreeGeneric *param) {
	char buffer[MAX_STRING_SIZE];
	fprintf((FILE*)param->info, "\\pswedge[linewidth=0.5pt,fillcolor=%s](%lf,%lf){%lfpt}{%lf}{%lf}\n", sprintRGBPSTricks(buffer, param->curgb), x, param->height-y, param->radius, 90*a/asin(1), 90*b/asin(1));
}

void fillGradientPSTricks(TypeRGB rgb0, TypeRGB rgb1, double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param) {
	char buffer1[MAX_STRING_SIZE], buffer2[MAX_STRING_SIZE];
	fprintf((FILE*)param->info,"\\psframe[linestyle=none,fillstyle=gradient,gradangle=90,gradbegin=%s,gradend=%s,gradmidpoint=1.](%.2lf,%.2lf)(%.2lf,%.2lf)",  sprintRGBPSTricks(buffer1, rgb0),  sprintRGBPSTricks(buffer2, rgb1), x0, param->height-y0, x1, param->height-y1);
}

void drawLineDotPSTricks(TypeRGB rgb, double alpha, double radius, double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param) {
	char buffer[MAX_STRING_SIZE];
	fprintf((FILE*)param->info, "\\psline[linewidth=%lf,linecolor=%s,opacity=%lf]{cc-cc}(%.2lf, %.2lf)(%.2lf, %.2lf)\n", radius, sprintRGBPSTricks(buffer, param->curgb), alpha, x0, y0, x1, y0);
}

void fillPolygonPSTricks(TypeRGB rgb, double alpha, double *x, double *y, int size, TypeParamDrawTreeGeneric *param) {
	char buffer[MAX_STRING_SIZE];
	int i;
	fprintf((FILE*)param->info, "\\pspolygon[fillstyle=solid,fillcolor=%s,opacity=%lf]", sprintRGBPSTricks(buffer, rgb), alpha);
	for(i=0; i<size; i++)
		fprintf((FILE*)param->info, "(%lf, %lf)", x[i], param->height-y[i]);
	fprintf((FILE*)param->info, "\n");
}

void drawDotPSTricks(TypeRGB rgb, double alpha, double radius, double x, double y, TypeParamDrawTreeGeneric *param) {
	char buffer[MAX_STRING_SIZE];
	fprintf((FILE*)param->info, "\\pscircle[linestyle=none,fillstyle=solid,fillcolor=%s,opacity=%lf](%.2lf, %.2lf){%.2lf}", sprintRGBPSTricks(buffer, rgb), alpha, x, y, radius);
}

double getMaxLeafLabelWidthPSTricks(TypeTree *tree) {
	int n;
	double max = 0.;
    if(tree->name)
		for(n=0; n<tree->size; n++)
            if(tree->node[n].child<0 && tree->name[n])
                if(strlen(tree->name[n])>max)
                    max = strlen(tree->name[n]);
	return max;
}

void setFunctPSTricks(TypeFunctDrawTreeGeneric *funct) {
	funct->start = startPSTricks;
	funct->end = endPSTricks;
	funct->drawText = drawTextPSTricks;
	funct->drawLine = drawLinePSTricks;
}

void startPSTricksStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param) {
	FILE *fo;
	if(!(fo = fopen(filename, "w"))) {
		fprintf(stderr, "Error while opening %s\n", filename);
		exit(1);
	}
	param->scale = 100.;
	param->xoffset = 3;
	param->yoffset = 0;
	param->leafSep = OFFSET;
	param->labelSep = 1.;
	param->ycenter = 3.;
	param->ydec = 2.;
	param->radius = 10.;
	param->roffset = 5.;
	param->leafCur = 0.;
	param->info = (void*) fo;
	param->height = height;
	param->width = width;
	param->tickLength = 3.;
    param->scale = (param->width-(param->xoffset+param->labelSep+param->labelWidth))/((param->tmax-param->tmin));
    param->curgb = (TypeRGB) {.red = 0., .green = 0., .blue = 0.};
	fprintf((FILE*)param->info, "\\begin{pspicture}(0,0)(%.2lf, %.2lf)\n", param->width, param->height);
}

void startPSTricks(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param) {
	FILE *fo;
	if(!(fo = fopen(filename, "w"))) {
		fprintf(stderr, "Error while opening %s\n", filename);
		exit(1);
	}
	param->scale = 100.;
	param->xoffset = 3;
	param->yoffset = 0;
	param->leafSep = OFFSET;
	param->labelSep = 1.;
	param->ycenter = 3.;
	param->ydec = 2.;
	param->radius = 10.;
	param->roffset = 5.;
	param->leafCur = 0.;
    param->tmin = tree->minTime;
    param->tmax = tree->maxTime;
	param->info = (void*) fo;
	param->height = param->leafSep*(countLeaves(tree)+1)+3*OFFSET+3*LABEL_SEP+FONT_SIZE;
	param->width = param->height*param->ratio;
	param->tickLength = 3.;
	param->labelWidth = getMaxLeafLabelWidthPSTricks(tree)*CHAR_WIDTH;
	if(param->tmin == param->tmax)
		param->tmax++;
    param->scale = (param->width-(param->xoffset+param->labelSep+param->labelWidth))/((param->tmax-param->tmin));
    param->curgb = (TypeRGB) {.red = 0., .green = 0., .blue = 0.};
	fprintf((FILE*)param->info, "\\begin{pspicture}(0,0)(%.2lf, %.2lf)\n", param->width, param->height);
}

void endPSTricks(TypeParamDrawTreeGeneric *param) {
	fprintf((FILE*)param->info, "\\end{pspicture}\n");
	fclose((FILE*)param->info);
}
