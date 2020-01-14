#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>

#include "Utils.h"
#include "Tree.h"
#include "Fossil.h"
#include "DrawTreeCairo.h"
#include "DrawTreePSTricks.h"
#include "DrawTreeTikz.h"
#include "DrawTreeGeneric.h"
#include "DrawFossilInt.h"


#ifdef DO_PS
#endif

#define STRING_SIZE 300
#define INC_CUT 5
#define SEQ_SIZE 30
#define EXT_OUTPUT "_added.phy"
#define MAX_PRED 7.5E+07

#define SIZE_BUFFER_CHAR 300
#define INFTY 1E99
#define RINFTY 1E99
#define DEF 10
#define MIN_VAL 0.000001
#define DELTA 0.000001

#define MINVAL 0.01
#define TRIAL 10
#define FIX_VAL(x) (((x)<=0)?MINVAL:(x))

#define MAX_ITER 1000

#define M_MAX 6
#define M_MAX_F 4
#define MIN_TREE 20
#define PREFIX "table"
#define MAX_TRIALS 1000
#define MAX_NODE 100

#define HELPMESSAGE "\n\nusage: Distribution [options] <input file> [<output file>]\n\nEstimate the diversification rates of the tree contained in the input file.\nThe input file has to be in Newick format with special tags for fossils ages and origin and end of the diversification, \nit returns a text report with the estimates.\n\nOptions are:\n\t-o <options file name>\tload the settings of the optimizer. <options file name> has to be in the format:\n\t\t:SPE [0;1] :EXT [0;1] :FOS [0:1] :TRI 10 :TOL 1.E-7 :ITE 500\n\t-h\tdisplay help\n\n"


int main(int argc, char **argv) {	
	char *inputFileNameTree, *inputFileNameFossil = NULL, *outputPrefix = PREFIX, outputDistribution[STRING_SIZE], option[256], format = '1';
	FILE *fi, *fo;
	int i, j;
	double minTimeIntInf = NO_TIME, minTimeIntSup = NO_TIME, maxTimeIntInf = 0., maxTimeIntSup = 0.;

	
	for(i=0; i<256; i++)
		option[i] = 0;
	   
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['x']) {
			option['x'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &format) == 1)
				i++;
			else
				exitProg(ErrorArgument, "a character is required after option -f");
		}
		if(option['z']) {
			option['z'] = 0;
			if((i+1)<argc && (fi = fopen(argv[++i], "r"))) {
				TypeTree *tree;
				tree = readTree(fi);
				toBinary(tree);
				printTreeDebug(stdout, tree->root, tree, tree->name);
//				bltoabsTime(tree);
				reorderTree(tree->name, tree);
				if(tree->minTime == NO_TIME || tree->minTime == 0.)
					tree->minTime = tree->time[tree->root]*0.9;
				if(tree->maxTime == NO_TIME) {
					int n;
					tree->maxTime = 0.;
					for(n=0; n<tree->size; n++)
						if(tree->time[n]>tree->maxTime)
							tree->maxTime = tree->time[n];
				}
				//TypeInfoDrawTreeGeneric info;
				//double *timeSave;
				//int n;
				//timeSave = tree->time;
				//tree->time = (double*) malloc(tree->size*sizeof(double));
				//for(n=0; n<tree->size; n++)
					//tree->time[n] = timeSave[n];
				//fillUnknownTimes(tree->minTime, tree->maxTime, tree);
				//setFunctPDF(&(info.funct));
				//drawTreeFileGenericDebug("tree.pdf", tree, &info, NULL);
				//free((void*) tree->time);
				//tree->time = timeSave;
			} else {
				fprintf(stderr, "Error while reading %s.\n", argv[i]);
				exit(1);
			}
			exit(0);
		}
		if(option['o']) {
			option['o'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &minTimeIntInf) == 1)
				i++;
			else
				exitProg(ErrorArgument, "at least one value is expected after -o");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &minTimeIntSup) == 1)
				i++;
			else
				minTimeIntSup = NO_TIME;
		}
		if(option['h']) {
			option['h'] = 0;
			printf("%s\n", HELPMESSAGE);
			exit(0);
		}
	}
	if(i<argc) {
		inputFileNameTree = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a phylogenetic tree in Newick format\n");
		exit(1);
	}
	if(i<argc) {
		inputFileNameFossil = argv[i++];
	} else {
		fprintf(stderr, "Please provide the name of a file containing a fossil list\n");
		exit(1);
	}
	if(i<argc)
		outputPrefix = argv[i++];
	if((fi = fopen(inputFileNameTree, "r"))) {
		TypeTree *tree;
		TypeFossilIntFeature *fos;
		TypeDistribution *d;
		int i, n;
		TypeInfoDrawTreeGeneric info;
		TypeAdditionalDrawTreeGeneric add;
		TypeDataDrawFossilInt data;
        tree = readTree(fi);
        fclose(fi);
		toBinary(tree);
		reorderTreeSize(tree);
		if(inputFileNameFossil != NULL) {
			FILE *ff;
			if((ff = fopen(inputFileNameFossil, "r"))) {
				fos = getFossilIntFeature(ff, tree->name, tree->size);
			} else {
				fprintf(stderr, "Cannot read %s\n", inputFileNameFossil);
				exit(1);
			}
		} else
			fos = fosToFossilInt(tree);
printf("Fossil %d %s\n", fos->sizeFossil, inputFileNameFossil);
		if(getMaxFossilIntTime(fos) > 0.)
			negateFossilInt(fos);
		fixStatus(tree, fos);
//fprintTreeFossilInt(stdout, tree, fos);
		if(nTodo>0)
			todo = tab_todo;
		else {
			int n;
			todo = (int*) malloc(sizeof(int)*tree->size/2);
			for(n=0; n<tree->size; n++)
				if(tree->node[n].child != NOSUCH)
					todo[nTodo++] = n;
		}
		for(n=0; n<tree->size; n++)
			if(tree->node[n].child != NOSUCH)
				tree->time[n] = NO_TIME;
		if(tree->minTimeInt.inf != NO_TIME)
			minTimeIntInf = tree->minTimeInt.inf;
		if(tree->minTimeInt.sup != NO_TIME)
			minTimeIntSup = tree->minTimeInt.sup;
		if(tree->maxTimeInt.inf != NO_TIME)
			maxTimeIntInf = tree->maxTimeInt.inf;
		if(tree->maxTimeInt.sup != NO_TIME)
			maxTimeIntSup = tree->maxTimeInt.sup;
		double minFossilTime = getMinFossilIntTime(fos);
		if(minTimeIntInf == NO_TIME && minTimeIntSup == NO_TIME) {
			double tmp = getMinFossilIntTime(fos);
			if(tmp<0) {
				minTimeIntInf = 1.2*minFossilTime;
				minTimeIntSup = 1.1*minFossilTime;
			} else {
				minTimeIntInf = 0.8*minFossilTime;
				minTimeIntSup = 0.9*minFossilTime;
			}
		}
		if(minTimeIntInf > minFossilTime)
			minTimeIntInf = minFossilTime;
		tree->maxTime = maxTimeIntSup;
		tree->minTime = minTimeIntInf;
		tree->minTimeInt.inf = minTimeIntInf;
		tree->minTimeInt.sup = minTimeIntSup;
		tree->maxTime = maxTimeIntSup;
		tree->maxTimeInt.inf = maxTimeIntInf;
		tree->maxTimeInt.sup = maxTimeIntSup;
		if(tree->parent==NULL)
			tree->parent = getParent(tree);
		}
		if(format != '0') {
			char *tmp, outputFileNameG[SIZE_BUFFER_CHAR], *outputFileName;
			outputFileName = inputFileNameTree;
			for(n=0; n<tree->size; n++) {
				if(tree->node[n].child == NOSUCH) {
							int f;
							double max;
					switch(fos->status[n]) {
						case contempNodeStatus:
							tree->time[n] = tree->maxTime;
						break;
						case unknownNodeStatus:
							//tree->time[n] = NO_TIME;
						break;
						case extinctNodeStatus:
						//case unknownNodeStatus:
							max = NEG_INFTY;
							for(f=fos->fossilInt[n]; f!=NOSUCH; f=fos->fossilIntList[f].prec)
								if(fos->fossilIntList[f].fossilInt.sup>max)
									max = fos->fossilIntList[f].fossilInt.sup; 
							tree->time[n] = max;
						break;
						default:
							fprintf(stderr, "Node %d has no status\n", n);
							return 1;
					}
				}
			}
			if((tmp = strrchr(outputFileName, '.')) != NULL)
				tmp[0] = '\0';
			info.param.tmin = tree->minTime;
			info.param.tmax = tree->maxTime;
			info.param.ratio = 0.55;
			info.param.scaleStep = 10.;
//			data.color = (TypeRGB) {.red = 1., .green = 1., .blue = 0.};
			data.color = (TypeRGB) {.red = 0.63, .green = 0.32, .blue = 0.18};
			data.radius = 5.;
			data.alpha = 0.5;
			data.fos = fos;
			dataD.color = (TypeRGB) {.red = 1., .green = 0., .blue = 0.};
			dataD.alpha = 0.5;
			dataD.dens = d;
			dataFD.dataFossil = &data;
			dataFD.dataDensity = &dataD;
			add.data = &data;
			add.draw = drawFossilInt;
			tree->maxTime = getMaximumTime(tree);
			fillUnknownTimesFossilInt(tree->minTime, tree->maxTime, tree, fos);
			switch(format) {
				case '1':
					sprintf(outputFileNameG, "%s_tree.pdf", outputFileName);
					setFunctPDF(&(info.funct));
					data.drawLineDot = drawLineDotCairo;
					dataD.fillPolygon = fillPolygonCairo;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '2':
					sprintf(outputFileNameG, "%s_tree.ps", outputFileName);
					setFunctPS(&(info.funct));
					data.drawLineDot = drawLineDotCairo;
					dataD.fillPolygon = fillPolygonCairo;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '3':
					sprintf(outputFileNameG, "%s_tree.png", outputFileName);
					setFunctPNG(&(info.funct));
					data.drawLineDot = drawLineDotCairo;
					dataD.fillPolygon = fillPolygonCairo;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '4':
					sprintf(outputFileNameG, "%s_tree.svg", outputFileName);
					setFunctSVG(&(info.funct));
					data.drawLineDot = drawLineDotCairo;
					dataD.fillPolygon = fillPolygonCairo;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '5':
					sprintf(outputFileNameG, "%s_tree_pst.tex", outputFileName);
					setFunctPSTricks(&(info.funct));
					data.drawLineDot = drawLineDotPSTricks;
					dataD.fillPolygon = fillPolygonPSTricks;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				case '6':
					sprintf(outputFileNameG, "%s_tree_tikz.tex", outputFileName);
					setFunctTikz(&(info.funct));
					data.drawLineDot = drawLineDotTikz;
					dataD.fillPolygon = fillPolygonTikz;
					drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
					break;
				default:
					;
			}
		}
		if(todo != NULL && todo != tab_todo)
			free((void*)todo);
		if(d != NULL) {
			for(i=0; i<tree->size; i++)
				if(d[i].item != NULL)
					free((void*)d[i].item);
			free((void*)d);
		}
		gsl_rng_free(rg);
		freeTree(tree);
		freeFossilIntFeature(fos);
	} else {
		fprintf(stderr, "Cannot read %s\n", inputFileNameTree);
		exit(1);
	}
	return 0;
}
