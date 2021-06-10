#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_rng.h>

#include "Utils.h"
#include "Tree.h"
#include "SimulTree.h"
#include "Fossil.h"
#include "SimulFossil.h"
#include "DrawTreeCairo.h"
#include "DrawTreePSTricks.h"
#include "DrawTreeTikz.h"
#include "DrawTreeGeneric.h"
#include "DrawFossil.h"
#include "DrawDensity.h"


#define STRING_SIZE 300
#define HELP_MESSAGE "\nusage: sample [options] [<output file>]\n\nsample simulates random trees and fossils finds and saves them in Newick format\n\nOptions:\n\t-h : display help\n\t-b <birth>\t: set birth rate\n\t-d <death>\t: set death rate\n\t-f <fossil>\t: set fossil find rate\n\t-m <min>\t: set minimum number of contemporary species of a simulation to be considered\n\t-M <size>\t: set maximum size of a simulation to be considered\n\t-i <niter>\t: set the number of simulations\n\t-t <time> : the end time of the diversification (start is always 0)\n"

//./samp -b 2 -d 1 -f 1 -t 5 -i 4 simul.newick
void nameLeavesRec(int n, int *ind, int length, TypeTree *tree);
void renameLeaves(TypeTree *tree);

int main(int argc, char **argv) {	
	char outputName[STRING_SIZE], outputFileName[STRING_SIZE], outputFileNameG[STRING_SIZE+50],*outputPrefix, option[256], format = '1';
	int i, nb, niter = 5, minContemp = 2, minFossil = 5, maxSizeTree = 10000, minSizeTree = 15;
	double birth = 2., death = 1., fossil = 1., maxTime = 5., figwidth = 100.;
	gsl_rng *rg = gsl_rng_alloc(gsl_rng_random_glibc2);	
	TypeAdditionalDrawTreeGeneric add;
	TypeDataDrawFossil data;

	for(i=0; i<256; i++)
		option[i] = 0;
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		int j;
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
		if(option['b']) {
			option['b'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &birth) == 1)
				i++;
		}
		if(option['d']) {
			option['d'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &death) == 1)
				i++;
		}
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &fossil) == 1)
				i++;
		}
		if(option['i']) {
			option['i'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &niter) == 1)
				i++;
		}
		if(option['m']) {
			option['m'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &minContemp) == 1)
				i++;
		}
		if(option['M']) {
			option['M'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &maxSizeTree) == 1)
				i++;
		}
		if(option['t']) {
			option['t'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &maxTime) == 1)
				i++;
		}
		if(option['x']) {
			option['x'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%c", &format) == 1)
				i++;
			else
				error("a character is required after option -f");
		}
		if(option['y']) {
			option['y'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &figwidth) == 1)
				i++;
			else
				error(ErrorArgument, "a number is required after option -t");
		}
		if(option['h']) {
			printf("%s\n", HELP_MESSAGE);
			exit(EXIT_SUCCESS);
		}
	}
	if(!(i<argc && sscanf(argv[i++], "%s", outputName) == 1)) {
		strcpy(outputName, "Tree_simul");
	}
	if((outputPrefix = strrchr(outputName, '.')) != NULL)
		outputPrefix[0] = '\0';
	if((outputPrefix=strrchr(outputName, '/')) == NULL)
		outputPrefix = outputName;
	else
		outputPrefix++;
	add.draw = drawFossil;
	data.color = (TypeRGB) {.red = 0.3, .green = 0., .blue = 0.};
	data.radius = 3.;
	data.alpha = 0.75;
	for(nb=1; nb <= niter; nb++) {
		TypeFossilFeature *fos;
		TypeTree *tree, *tree1, *tree2;
		FILE *ft;
		char outputName[STRING_SIZE];
		TypeInfoDrawTreeGeneric info;
		tree = NULL;
		do {
			if(tree != NULL)
				freeTree(tree);
			tree = simulTree(rg, birth, death, maxTime);
//		} while(!(tree != NULL && (countContemp(tree))>=minContemp && tree->size>=minSizeTree && tree->size<=maxSizeTree));
		} while(!(tree != NULL && (countContemp(tree))==0 && tree->size>=minSizeTree && tree->size<=maxSizeTree));
		reorderTreeSize(tree);
		renameLeaves(tree);
		fos = NULL;
		do {
			if(fos != NULL)
				freeFossilFeature(fos);
			fos = addFossils(rg, fossil, tree);
		} while(fos->size <= minFossil);
		info.param.tmin = tree->minTime;
		info.param.tmax = tree->maxTime;
		info.param.scaleStep = 0.1;
		info.param.width = figwidth;
		data.fos = fos;
		add.data = (void*) &data;
		sprintf(outputFileName, "%s_sim_all_%d", outputPrefix, nb);
		switch(format) {
			case '1':
				sprintf(outputFileNameG, "%s_tree.pdf", outputFileName);
				setFunctPDF(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, NULL);
				break;
			case '2':
				sprintf(outputFileNameG, "%s_tree.ps", outputFileName);
				setFunctPS(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, NULL);
				break;
			case '3':
				sprintf(outputFileNameG, "%s_tree.png", outputFileName);
				setFunctPNG(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, NULL);
				break;
			case '4':
				sprintf(outputFileNameG, "%s_tree.svg", outputFileName);
				setFunctSVG(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, NULL);
				break;
			case '5':
				sprintf(outputFileNameG, "%s_tree_pst.tex", outputFileName);
				setFunctPSTricks(&(info.funct));
				data.drawDot = drawDotPSTricks;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, NULL);
				break;
			case '6':
				sprintf(outputFileNameG, "%s_tree_tikz.tex", outputFileName);
				setFunctTikz(&(info.funct));
				data.drawDot = drawDotTikz;
				data.radius = 0.08;
				drawTreeFileGeneric(outputFileNameG, tree, &info, &add, NULL);
				break;
			default:
				;
		}
		tree->name = NULL;
		tree1 = pruneFossilBis(tree, fos);
		freeFossilFeature(fos);
		freeTree(tree);
		tree2 = fixBinaryFossil(tree1, (TypeFossilFeature*) tree1->info);
		freeFossilFeature((TypeFossilFeature*) tree1->info);
		freeTree(tree1);
		sprintf(outputFileName, "%s_sim_obs_%d", outputPrefix, nb);
		data.fos = (TypeFossilFeature*) tree2->info;
		switch(format) {
			case '1':
				sprintf(outputFileNameG, "%s_tree.pdf", outputFileName);
				setFunctPDF(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '2':
				sprintf(outputFileNameG, "%s_tree.ps", outputFileName);
				setFunctPS(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '3':
				sprintf(outputFileNameG, "%s_tree.png", outputFileName);
				setFunctPNG(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '4':
				sprintf(outputFileNameG, "%s_tree.svg", outputFileName);
				setFunctSVG(&(info.funct));
				data.drawDot = drawDotCairo;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '5':
				sprintf(outputFileNameG, "%s_tree_pst.tex", outputFileName);
				setFunctPSTricks(&(info.funct));
				data.drawDot = drawDotPSTricks;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			case '6':
				sprintf(outputFileNameG, "%s_tree_tikz.tex", outputFileName);
				setFunctTikz(&(info.funct));
				data.drawDot = drawDotTikz;
				data.radius = 0.08;
				drawTreeFileGeneric(outputFileNameG, tree2, &info, &add, NULL);
				break;
			default:
				;
		}
		sprintf(outputName, "%s_sim_%d.newick", outputPrefix, nb);
		if((ft = fopen(outputName, "w"))) {
			int n;
			fprintf(ft, "\n\n");
			tree2->minTimeInt.inf = 0.;
			tree2->minTimeInt.sup = 0.;
			tree2->maxTimeInt.inf = maxTime;
			tree2->maxTimeInt.sup = maxTime;
			tree2->minTime = 0.;
			tree2->maxTime = maxTime;
			tree2->minTimeInt.inf = 0.;
			tree2->minTimeInt.sup = 0.;
			tree2->maxTimeInt.inf = maxTime;
			tree2->maxTimeInt.sup = maxTime;
			//for(n=0; n<tree2->size; n++)
				//if(tree2->time[n] != maxTime)
					//tree2->time[n] = NO_TIME;
			tree2->name = nameLeaves("leaf_", tree2);
			tree2->comment = (char**) malloc(tree2->sizeBuf*sizeof(char*));
			for(n=0; n<tree2->sizeBuf; n++)
				tree2->comment[n] = NULL;
			fillCommentFossil(tree2, (TypeFossilFeature*) tree2->info);
			fprintSubtreeNewick(ft, tree2->root, tree2);
			fprintf(ft, "\n");
			fclose(ft);
		}
		freeFossilFeature((TypeFossilFeature*) tree2->info);
		freeTree(tree2);
	}
	return 0;
}

void nameLeavesRec(int n, int *ind, int length, TypeTree *tree) {
	if(tree->node[n].child < 0) {
		int i, tmp;
		tree->name[n] = (char*) malloc((length+1)*sizeof(char));
		tmp = *ind;
		tree->name[n][length] = '\0';
		for(i=0; i<length; i++) {
			tree->name[n][length-1-i] = '0'+tmp%10;
			tmp /= 10;
		}
		(*ind)++;
	}
}

void renameLeaves(TypeTree *tree) {
	int length, ind;
	if(tree->name == NULL) {
		int i;
		tree->name = (char**) malloc(tree->sizeBuf*sizeof(char*));
		for(i=0; i<tree->sizeBuf; i++)
			tree->name[i] = NULL;
	}
	length = (int) ceil(log(tree->size)/log(10.));
	ind = 0;
	nameLeavesRec(tree->root, &ind, length, tree);
}
