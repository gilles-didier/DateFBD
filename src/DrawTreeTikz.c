#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <locale.h>
#include "Utils.h"
#include "DrawTreeTikz.h"

#define NB_MIN 6
#define OFFSET 0.25
#define TICK_LENGTH 0.1
#define FONT_NAME "Helvetica"
//#define LABEL_SEP 10.
#define MAX_STRING_SIZE 200
//#define STANDARD_WIDTH 500
#define CHAR_WIDTH 0.1
#define LABEL_SEP 0.2
#define FONT_SIZE 9.
#define STANDARD_WIDTH 500
#define MY_PI acos(-1)



char *sprintRGBTikz(char *buffer, TypeRGB rgb) {
	sprintf(buffer, "\\definecolor{color}{rgb}{%.3lf,%.3lf,%.3lf}", rgb.red, rgb.green, rgb.blue);
	return buffer;
}

void drawTextTikz(double x0, double y0, char *text, char *mod, TypeParamDrawTreeGeneric *param) {
    int i;
    char *tmp = strdpl(text), or[50];
    for(i=0; tmp[i]!='\0'; i++)
        if(tmp[i] == '_')
            tmp[i] = ' ';
    if(mod == NULL || mod[0] == '\0')
		strcpy(or, "center");
	else {
		switch(mod[0]) {
			case 't':
				strcpy(or, "north");
				break;
			case 'l':
				strcpy(or, "west");
				break;
			case 'r':
				strcpy(or, "east");
				break;
			case 'b':
				strcpy(or, "south");
				break;
			default:
				strcpy(or, "center");
		}
		switch(mod[1]) {
			case 't':
				strcat(or, "north");
				break;
			case 'l':
				strcat(or, "west");
				break;
			case 'r':
				strcat(or, "east");
				break;
			case 'b':
				strcat(or, "south");
				break;
			default:
				;
		}
	}
    fprintf((FILE*)param->info, "\\node[anchor=%s] at (%.2lf, %.2lf) {\\small \\textcolor[rgb]{%.3lf %.3lf %.3lf}{\\textit{%s}}};\n", or, x0, param->height-y0, param->curgb.red, param->curgb.green, param->curgb.blue, tmp);
    free((void*)tmp);
}

void drawTextAngleTikz(double x0, double y0, double a, char *text, char *mod, TypeParamDrawTreeGeneric *param) {
    int i;
    char *tmp = strdpl(text), *or;
    for(i=0; tmp[i]!='\0'; i++)
        if(tmp[i] == '_')
            tmp[i] = ' ';
    if(mod == NULL || mod[0] == '\0')
		or = "center";
	else
		switch(mod[0]) {
			case 't':
				or = "north";
				break;
			case 'l':
				or = "west";
				break;
			case 'r':
				or = "east";
				break;
			case 'b':
				or = "south";
				break;
			default:
				or = "center";
		}
    fprintf((FILE*)param->info, "\\node[rotate=%.2lf, anchor=%s] at (%.2lf, %.2lf) {\\small \\textcolor[rgb]{%.3lf %.3lf %.3lf}{%s}};\n", -(180.*a)/MY_PI, or, x0, param->height-y0, param->curgb.red, param->curgb.green, param->curgb.blue, tmp);
//	fprintf((FILE*)param->info, "\\rput[%s]{%.2lf}(%.2lf,%.2lf){\\small \\textcolor[rgb]{%.3lf %.3lf %.3lf}{%s}}\n", mod, -(180.*a)/MY_PI, x0, param->height-y0, param->curgb.red, param->curgb.green, param->curgb.blue, tmp);
	free((void*)tmp);
}

void drawLineTikz(double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param) {
	char buffer[MAX_STRING_SIZE];
    fprintf((FILE*)param->info, "%s\n\\draw[line width=1pt,draw=color] (%.2lf, %.2lf) -- (%.2lf, %.2lf);\n", sprintRGBTikz(buffer, param->curgb), x0, param->height-y0, x1, param->height-y1);
//    fprintf((FILE*)param->info, "\\psline[line width=1pt,draw=%s](%.2lf, %.2lf)(%.2lf, %.2lf)\n", sprintRGBTikz(buffer, param->curgb), x0, param->height-y0, x1, param->height-y1);
}

void drawLineDotTikz(TypeRGB rgb, double alpha, double radius, double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param) {
	char buffer[MAX_STRING_SIZE];
    fprintf((FILE*)param->info, "%s\n\\draw[line width=%.2lf,opacity=%.2lf,draw=color,line cap=round] (%.2lf, %.2lf) -- (%.2lf, %.2lf);\n", sprintRGBTikz(buffer, rgb), radius, alpha, x0, param->height-y0, x1, param->height-y1);
//    fprintf((FILE*)param->info, "\\psline[line width=1pt,draw=%s](%.2lf, %.2lf)(%.2lf, %.2lf)\n", sprintRGBTikz(buffer, param->curgb), x0, param->height-y0, x1, param->height-y1);
}

void drawDottedLineTikz(double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param) {
//	char buffer[MAX_STRING_SIZE];
//    fprintf((FILE*)param->info, "\\psline[linestyle=dotted,line width=1pt,draw=%s](%.2lf, %.2lf)(%.2lf, %.2lf)\n", sprintRGBTikz(buffer, param->curgb), x0, param->height-y0, x1, param->height-y1);
}

void fillWedgeTikz(TypeRGB rgb, double x, double y, double a, double b, TypeParamDrawTreeGeneric *param) {
//	char buffer[MAX_STRING_SIZE];
//	fprintf((FILE*)param->info, "\\pswedge[linestyle=none,fillstyle=solid,fillcolor=%s,fillcolor=%s](%lf,%lf){%lfpt}{%lf}{%lf}\n", sprintRGBTikz(buffer, rgb), sprintRGBTikz(buffer, param->curgb), x, param->height-y, param->radius, 90*a/asin(1), 90*b/asin(1));
}

void drawWedgeTikz(double x, double y, double a, double b, TypeParamDrawTreeGeneric *param) {
//	char buffer[MAX_STRING_SIZE];
//	fprintf((FILE*)param->info, "\\pswedge[line width=0.5pt,fillcolor=%s](%lf,%lf){%lfpt}{%lf}{%lf}\n", sprintRGBTikz(buffer, param->curgb), x, param->height-y, param->radius, 90*a/asin(1), 90*b/asin(1));
}

void fillGradientTikz(TypeRGB rgb0, TypeRGB rgb1, double x0, double y0, double x1, double y1, TypeParamDrawTreeGeneric *param) {
//	char buffer1[MAX_STRING_SIZE], buffer2[MAX_STRING_SIZE];
//	fprintf((FILE*)param->info,"\\psframe[linestyle=none,fillstyle=gradient,gradangle=90,gradbegin=%s,gradend=%s,gradmidpoint=1.](%.2lf,%.2lf)(%.2lf,%.2lf)",  sprintRGBTikz(buffer1, param->start),  sprintRGBTikz(buffer2, param->end), x0, param->height-y0, x1, param->height-y1);
}

void fillPolygonTikz(TypeRGB rgb, double alpha, double *x, double *y, int size, TypeParamDrawTreeGeneric *param) {
	int i;
	char buffer[MAX_STRING_SIZE];
	fprintf((FILE*)param->info, "%s\n\\fill[fill=color, fill opacity=%lf]", sprintRGBTikz(buffer, rgb), alpha);
	for(i=0; i<size; i++)
		fprintf((FILE*)param->info, "(%lf, %lf) -- ", x[i], param->height-y[i]);
	fprintf((FILE*)param->info, "cycle;\n");
}

void drawDotTikz(TypeRGB rgb, double alpha, double radius, double x, double y, TypeParamDrawTreeGeneric *param) {
	char buffer[MAX_STRING_SIZE];
	fprintf((FILE*)param->info, "%s\n\\fill[fill=color, fill opacity=%lf] (%lf, %lf) circle (%lfcm);\n", sprintRGBTikz(buffer, rgb), alpha, x, param->height-y, radius);
}

double getMaxLeafLabelWidthTikz(TypeTree *tree) {
	int n;
	double max = 0.;
    if(tree->name)
		for(n=0; n<tree->size; n++)
            if(tree->node[n].child<0 && tree->name[n])
                if(strlen(tree->name[n])>max)
                    max = strlen(tree->name[n]);
	return max;
}

void setFunctTikz(TypeFunctDrawTreeGeneric *funct) {
	funct->start = startTikz;
	funct->end = endTikz;
	funct->drawText = drawTextTikz;
	funct->drawLine = drawLineTikz;
}

void startTikzStd(char *filename, double width, double height, TypeParamDrawTreeGeneric *param) {
	FILE *fo;
	if(!(fo = fopen(filename, "w"))) {
		fprintf(stderr, "Error while opening %s\n", filename);
		exit(1);
	}
	param->scale = 10.;
	param->xoffset = 0.1;
	param->yoffset = 0;
	param->leafSep = OFFSET;
	param->labelSep = 0.1;
	param->ycenter = 0.1;
	param->ydec = 0.1;
	param->radius = 5.;
	param->roffset = 5.;
	param->leafCur = 0.;
	param->info = (void*) fo;
	param->height = height;
	param->width = height*0.75;
	param->tickLength = 0.02;
    param->scale = (param->width-(param->xoffset+param->labelSep+param->labelWidth))/((param->tmax-param->tmin));
    param->curgb = (TypeRGB) {.red = 0., .green = 0., .blue = 0.};
	fprintf((FILE*)param->info, "\\begin{tikzpicture}\n");
}

void startTikz(char *filename, TypeTree *tree, TypeParamDrawTreeGeneric *param) {
	FILE *fo;
	if(!(fo = fopen(filename, "w"))) {
		fprintf(stderr, "Error while opening %s\n", filename);
		exit(1);
	}
	param->scale = 10.;
	param->xoffset = 0.1;
	param->yoffset = 0;
	param->leafSep = OFFSET;
	param->labelSep = 0.01;
	param->ycenter = 0.05;
	param->ydec = 0.05;
	param->radius = 5.;
	param->roffset = 0.05;
	param->leafCur = 0.;
    param->tmin = tree->minTime;
    param->tmax = tree->maxTime;
	param->info = (void*) fo;
	param->height = param->leafSep*(countLeaves(tree)+1)+3*OFFSET+3*LABEL_SEP+FONT_SIZE;
	param->width = param->height*param->ratio;
	param->tickLength = 0.07;
	param->labelWidth = getMaxLeafLabelWidthTikz(tree)*CHAR_WIDTH;
	if(param->tmin == param->tmax)
		param->tmax++;
    param->scale = (param->width-(param->xoffset+param->labelSep+param->labelWidth))/((param->tmax-param->tmin));
    param->curgb = (TypeRGB) {.red = 0., .green = 0., .blue = 0.};
	fprintf((FILE*)param->info, "\\begin{tikzpicture}\n");
}

void endTikz(TypeParamDrawTreeGeneric *param) {
	fprintf((FILE*)param->info, "\\end{tikzpicture}\n");
	fclose((FILE*)param->info);
}
