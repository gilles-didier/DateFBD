#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <locale.h>
#include "Utils.h"
#include "DrawTreeGeneric.h"

#define NB_MIN 6
#define OFFSET 10
#define TICK_LENGTH 10
#define FONT_NAME "Helvetica"
#define FONT_SIZE 9.
#define LABEL_SEP 10.
#define MAX_STRING_SIZE 500
#define STANDARD_WIDTH 1000
#define CHAR_WIDTH 10

static double drawNodeGeneric(int n, int parent, TypeTree *tree,  TypeInfoDrawTreeGeneric *param, TypeAdditionalDrawTreeGeneric *add, double *ytab);
static void drawTreeGeneric(TypeTree *tree, TypeInfoDrawTreeGeneric *param, TypeAdditionalDrawTreeGeneric *add, TypeAdditionalDrawTreeGenericGeneral *addG);


void drawScaleGeneric(double x, double y, TypeInfoDrawTreeGeneric *info) {
	int flag = 0;
	double start, end, step, cur, width;
	width = info->param.width-info->param.labelWidth-info->param.xoffset-info->param.labelSep;
printf("line %.2lf %.2lf %.2lf\n", y, x, x+width);
	info->funct.drawLine(x, y, x+width, y, &(info->param));
	if((info->param.tmax-info->param.tmin) <= 0.)
		return;
	step = info->param.scaleStep;
	//step = pow(10., floor(log10(info->param.tmax-info->param.tmin)));
	//if((info->param.tmax-info->param.tmin)/step<NB_MIN)
		//step /= 2.;
	flag = step<1.;
	start = step*ceil(info->param.tmin/step);
	end = step*floor(info->param.tmax/step);
//printf("step %.1lf, max %.1lf\n", step, info->param.tmax);
	for(cur=start; cur<=end; cur += step) {
		char tmp[500];
		info->funct.drawLine(x+(cur-info->param.tmin)*info->param.scale, y, x+(cur-info->param.tmin)*info->param.scale, y+info->param.tickLength, &(info->param));
		if(flag)
			sprintf(tmp, "%.1lf", cur);
		else
			sprintf(tmp, "%.0lf", cur);
		info->funct.drawText(x+(cur-info->param.tmin)*info->param.scale, y+info->param.tickLength, tmp, "t", &(info->param));
	}
}

void drawTreeFileGeneric(char *filename, TypeTree *tree, TypeInfoDrawTreeGeneric *info, TypeAdditionalDrawTreeGeneric *add, TypeAdditionalDrawTreeGenericGeneral *addG) {
	int i;
	double *timeSave;
	
	timeSave = tree->time;
	tree->time = (double*) malloc(tree->size*sizeof(double));
	if(tree->time != NULL) {
		int hasUnknown = 0;
		for(i=0; i<tree->size; i++) {
			tree->time[i] = timeSave[i];
			if(timeSave[i] == NO_TIME)
				hasUnknown = 1;
		}
		if(hasUnknown)
			fillUnknownTimes(info->param.tmin, info->param.tmax,  tree);
	} else {
		for(i=0; i<tree->size; i++)
			tree->time[i] = 1.;
		bltoabsTime(tree);
	}
//	info->param.nleaves = countLeaves(tree);
	info->param.nleaves = tree->size/2+1;
printf("ok %s\n", filename);
	info->funct.start(filename, tree, &(info->param));
printf("tree->height %.2lf\n", info->param.height);
printf("tree scale %.2lf\n", info->param.leafSep*(countLeaves(tree)+1)+info->param.labelSep);
printf(" %d\n", tree->size);
printf("ok\n");
	drawTreeGeneric(tree, info, add, addG);
printf("ok\n");
	drawScaleGeneric(info->param.xoffset, info->param.leafSep*(info->param.nleaves+1)+info->param.labelSep, info);
printf("ok\n");
	info->funct.end(&(info->param));
	free((void*) tree->time);
	tree->time = timeSave;
}

void drawTreeFileGenericDebug(char *filename, TypeTree *tree, TypeInfoDrawTreeGeneric *info, TypeAdditionalDrawTreeGeneric *add) {
	int i;
	double *timeSave;
	char **nameSave;
	
	timeSave = tree->time;
	tree->time = (double*) malloc(tree->size*sizeof(double));
	if(tree->time != NULL) {
		int hasUnknown = 0;
		for(i=0; i<tree->size; i++) {
			tree->time[i] = timeSave[i];
			if(timeSave[i] == NO_TIME)
				hasUnknown = 1;
		}
		if(hasUnknown)
			fillUnknownTimes(info->param.tmin, info->param.tmax,  tree);
	} else {
		for(i=0; i<tree->size; i++)
			tree->time[i] = 1.;
		bltoabsTime(tree);
	}
	nameSave = tree->name;
	tree->name = (char**) malloc(tree->size*sizeof(char*));
	for(i=0; i<tree->size; i++) {
		char buffer[200];
		if(nameSave != NULL && nameSave[i] != NULL)
			sprintf(buffer, "%s/%d", nameSave[i], i);
		else
			sprintf(buffer, "%d", i);
		tree->name[i] = (char*) malloc((strlen(buffer)+1)*sizeof(char));
		strcpy(tree->name[i], buffer);
	}
	info->funct.start(filename, tree, &(info->param));
	drawTreeGeneric(tree, info, add, NULL);
	drawScaleGeneric(info->param.xoffset, info->param.leafSep*(countLeaves(tree)+1)+info->param.labelSep, info);
	info->funct.end(&(info->param));
	for(i=0; i<tree->size; i++)
		if(tree->name[i] != NULL)
			free((void*)tree->name[i]);
	free((void*) tree->name);
	tree->name = nameSave;
	free((void*) tree->time);
	tree->time = timeSave;
}

void drawTreeGeneric(TypeTree *tree, TypeInfoDrawTreeGeneric *info, TypeAdditionalDrawTreeGeneric *add, TypeAdditionalDrawTreeGenericGeneral *addG) {
	int tmp, n;
	double min, max, *ytab;
	if(tree->size<=0)
		return;
	ytab = (double*) malloc(tree->size*sizeof(double));
	if((tmp = tree->node[tree->root].child) >= 0) {
		min = drawNodeGeneric(tmp, tree->root, tree, info, add, ytab);
		max = min;
		for(tmp = tree->node[tmp].sibling; tmp >= 0; tmp = tree->node[tmp].sibling)
			max = drawNodeGeneric(tmp, tree->root, tree, info, add, ytab);
	} else {
		max = info->param.leafCur+info->param.leafSep/2.;
		min = max;
	}
	info->funct.drawLine((tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset, min, (tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset, max, &(info->param));
	ytab[tree->root] = (min+max)/2;
	info->funct.drawLine((tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+ytab[tree->root], info->param.xoffset, info->param.yoffset+ytab[tree->root], &(info->param));
	if(tree->name != NULL && tree->name[tree->root] != NULL)
		info->funct.drawText((tree->time[tree->root]-info->param.tmin)*info->param.scale+info->param.xoffset+info->param.labelSep, ytab[tree->root]+info->param.yoffset, tree->name[tree->root], "l", &(info->param));
	if(add != NULL)
		for(n=0; n<tree->size; n++)
			add->draw(n, (tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, ytab[n], info, add->data);
	if(addG != NULL)
		addG->draw(info, addG->data);
	free((void*)ytab);
}

double drawNodeGeneric(int n, int parent, TypeTree *tree, TypeInfoDrawTreeGeneric *info, TypeAdditionalDrawTreeGeneric *add, double *ytab) {
	double min, max, y;
	if(tree->node[n].child != NOSUCH) {
		int tmp = tree->node[n].child;
		min = drawNodeGeneric(tmp, n, tree, info, add, ytab);
		max = min;
		for(tmp = tree->node[tmp].sibling; tmp != NOSUCH; tmp = tree->node[tmp].sibling)
			max = drawNodeGeneric(tmp, n, tree, info, add, ytab);
		y = (min+max)/2;
		info->funct.drawLine((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+min, (tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+max, &(info->param));
		if(tree->name && tree->name[n])
			info->funct.drawText((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset+info->param.labelSep, y+info->param.yoffset, tree->name[n], "l", &(info->param));
	} else {
		info->param.leafCur += info->param.leafSep;
		y = info->param.leafCur;
		if(tree->name && tree->name[n])
			info->funct.drawText((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset+info->param.labelSep, info->param.yoffset+y, tree->name[n], "l", &(info->param));
	}
	info->funct.drawLine((tree->time[n]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, (tree->time[parent]-info->param.tmin)*info->param.scale+info->param.xoffset, info->param.yoffset+y, &(info->param));
	ytab[n] = y;
	return y;
}
















