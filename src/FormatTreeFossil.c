#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "Utils.h"
#include "Tree.h"
#include "Fossil.h"
#include "FossilInt.h"




#define HELPMESSAGE "\n\nusage: estimate [options] <input file> [<output file>]\n\nEstimate the diversification rates of the tree contained in the input file.\nThe input file has to be in Newick format with special tags for fossils ages and origin and end of the diversification, \nit returns a text report with the estimates.\n\nOptions are:\n\t-o <options file name>\tload the settings of the optimizer. <options file name> has to be in the format:\n\t\t:SPE [0;1] :EXT [0;1] :FOS [0:1] :TRI 10 :TOL 1.E-7 :ITE 500\n\t-h\tdisplay help\n\n"



int main(int argc, char **argv) {	
	char *inputFileNameTree, *inputFileNameFossil, *outputFileName, option[256];
	FILE *fi, *fo, *ff;
	int i, j;
	double minTimeIntInf = NO_TIME, minTimeIntSup = NO_TIME, maxTimeIntInf = 0., maxTimeIntSup = 0.;

	for(i=0; i<256; i++)
	option[i] = 0;

	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
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
		outputFileName = argv[i++];
	else
		outputFileName = "out.phy";
	if((fi = fopen(inputFileNameTree, "r")) && (ff = fopen(inputFileNameFossil, "r"))) {
		TypeTree **tree;
		TypeFossilIntFeature **fos;
		int sizeTree;
		tree = readTrees(fi);
		fclose(fi);
		if(tree[0] == NULL) {
			fprintf(stderr, "Error: no tree\n");
			exit(1);
		}
		for(sizeTree=0; tree[sizeTree]!=NULL; sizeTree++) {
			if(tree[sizeTree]->name!=NULL)
				for(i=0; i<tree[sizeTree]->size; i++)
					if(tree[sizeTree]->name[i]!=NULL)
						fixSpace(tree[sizeTree]->name[i]);
		}	
		fos = (TypeFossilIntFeature**) malloc(sizeTree*sizeof(TypeFossilIntFeature*));
		if((fo = fopen(outputFileName, "w"))) {
			for(i=0; i<sizeTree; i++) {
				double minFossilTime;
				toBinary(tree[i]);
				rewind(ff);
				fos[i] = getFossilIntFeature(ff, tree[i]->name, tree[i]->size);
				if(getMaxFossilIntTime(fos[i]) > 0.) {
					negateFossilInt(fos[i]);
				}
				fixStatus(tree[i], fos[i]);
				minFossilTime = getMinFossilIntTime(fos[i]);
				if(minTimeIntInf == NO_TIME || minTimeIntInf > minFossilTime) {
					if(minFossilTime<0)
						minTimeIntInf = 1.2*minFossilTime;
					else
						minTimeIntInf = 0.8*minFossilTime;
				}
				tree[i]->minTime = minTimeIntInf;
				tree[i]->minTimeInt.inf = minTimeIntInf;
				tree[i]->minTimeInt.sup = minTimeIntSup;
				tree[i]->maxTime = maxTimeIntInf;
				tree[i]->maxTimeInt.inf = maxTimeIntInf;
				tree[i]->maxTimeInt.sup = maxTimeIntSup;
				fprintTreeFossilInt(fo, tree[i], fos[i]);
				fprintf(fo, "\n\n");
			}
			fclose(fo);
		}
	} else {
		fprintf(stderr, "Cannot read %s or %s\n", inputFileNameTree, inputFileNameFossil);
		exit(1);
	}
	return 0;
}
