#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "Utils.h"
#include "Tree.h"
#include "Fossil.h"
#include "FossilInt.h"
#include "DrawTreeCairo.h"
#include "DrawTreePSTricks.h"
#include "DrawTreeTikz.h"
#include "DrawTreeGeneric.h"
#include "DrawFossilInt.h"
#include "DrawDensity.h"

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
#define PREFIX "outTree"
#define MAX_TRIALS 1000
#define MAX_NODE 100

#define HELPMESSAGE "--------------------------\n\nNAME\n	draw - Drawing a single tree with the fossil ages\n	\nSYNOPSIS\n	draw [OPTIONS] <input Tree File> <input Fossil File>\n\nDESCRIPTION\n	Output a figure in various graphic formats of the tree in <input Tree File> with the fossil ages of  <input Fossil File>\n\n	Options are\n	-z <input Tree File>\n		output the tree in text debug format in the console and exit \n	-o <origin> \n		set the origin time\n	-e <end> \n		set the end time\n	-x <number>\n		set the graphic format of the output (option is required if one wants a graphic output)\n			-f 1 -> pdf\n			-f 2 -> postscript\n			-f 3 -> png\n			-f 4 -> svg\n			-f 5 -> LaTeX (psTricks)\n			-f 6 -> LaTeX (TikZ)\n	-h\n		display help\n\n--------------------------\n"




int main(int argc, char **argv) {	
	char *inputFileNameTree, *inputFileNameFossil = NULL, *tmp, outputFileNameG[SIZE_BUFFER_CHAR], *outputFileName, option[256], format = '1';
	FILE *fi;
	int i, j;
	double minTimeIntInf = NO_TIME, maxTimeIntSup = 0.;

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
				printTreeDebug(stdout, tree->root, tree, tree->name);
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
				exitProg(ErrorArgument, "1 values are expected after -o");
		}
		if(option['e']) {
			option['e'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &maxTimeIntSup) == 1)
				i++;
			else
				exitProg(ErrorArgument, "1 values are expected after -o");
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
		fprintf(stderr, "Warning: No name of file containing a fossil list was provided\n");
	}
	if((fi = fopen(inputFileNameTree, "r"))) {
		TypeTree *tree;
		TypeFossilIntFeature *fos;
		int n;
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
		if(getMaxFossilIntTime(fos) > 0.)
			negateFossilInt(fos);
		fixStatus(tree, fos);
		double minFossilTime = getMinFossilIntTime(fos);
		if(minTimeIntInf == NO_TIME) {
			double tmp = getMinFossilIntTime(fos);
			if(tmp<0) {
				minTimeIntInf = 1.2*minFossilTime;
			} else {
				minTimeIntInf = 0.8*minFossilTime;
			}
		}
		if(minTimeIntInf > minFossilTime)
			minTimeIntInf = minFossilTime;
		tree->maxTime = maxTimeIntSup;
		tree->minTime = minTimeIntInf;
		tree->minTimeInt.inf = minTimeIntInf;
		tree->minTimeInt.sup = minTimeIntInf;
		if(tree->parent==NULL)
			tree->parent = getParent(tree);
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
		info.param.tmax = getMaximumTime(tree);
		info.param.ratio = 0.71;
		info.param.scaleStep = 10.;
		data.color = (TypeRGB) {.red = 0.63, .green = 0.32, .blue = 0.18};
		data.radius = 2.;
		data.alpha = 0.5;
		data.fos = fos;
		add.data = &data;
		add.draw = drawFossilInt;
		tree->maxTime = info.param.tmax;
		switch(format) {
			case '1':
				sprintf(outputFileNameG, "%s_tree.pdf", outputFileName);
				setFunctPDF(&(info.funct));
				data.drawLineDot = drawLineDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
				break;
			case '2':
				sprintf(outputFileNameG, "%s_tree.ps", outputFileName);
				setFunctPS(&(info.funct));
				data.drawLineDot = drawLineDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
				break;
			case '3':
				sprintf(outputFileNameG, "%s_tree.png", outputFileName);
				setFunctPNG(&(info.funct));
				data.drawLineDot = drawLineDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
				break;
			case '4':
				sprintf(outputFileNameG, "%s_tree.svg", outputFileName);
				setFunctSVG(&(info.funct));
				data.drawLineDot = drawLineDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
				break;
			case '5':
				sprintf(outputFileNameG, "%s_tree_pst.tex", outputFileName);
				setFunctPSTricks(&(info.funct));
				data.drawLineDot = drawLineDotPSTricks;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
				break;
			case '6':
				sprintf(outputFileNameG, "%s_tree_tikz.tex", outputFileName);
				setFunctTikz(&(info.funct));
				data.drawLineDot = drawLineDotTikz;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add);
				break;
			default:
				;
		}
		freeTree(tree);
		freeFossilIntFeature(fos);
	} else {
		fprintf(stderr, "Cannot read %s or %s\n", inputFileNameTree, inputFileNameFossil);
		exit(1);
	}
	return 0;
}
