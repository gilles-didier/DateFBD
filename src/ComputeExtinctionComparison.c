#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <fenv.h>


#include <gsl/gsl_rng.h>

#include "Utils.h"
#include "Tree.h"
#include "Fossil.h"
#include "FossilInt.h"
#include "Model.h"
#include "Uncertainty.h"
#include "MCMCImportanceSampling.h"

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
#define MAX_NODE 500
#define MAX_CLADE 500

#define HELPMESSAGE "--------------------------\n\nNAME\n	dist - Computation of the divergence time distibution associated to a given clade from a set of possible trees, the fossil ages, and the diversification rates\n	\nSYNOPSIS\n	dist [OPTIONS] <input Tree File> <input Fossil File> <input List Clade> [output File]\n\nDESCRIPTION\n	Compute the distribution of the divergence time associated to the clade corresponding to the list of tips given in <input List Clade> by sampling into the trees contianed in <input Tree File> (in Newick format) into the fossil ranges provided in <input Fossil File> (in csv format) and output the distribution as a .csv table <output File>.csv\n\n	Options are\n	-z <input Tree File>\n		output the tree in text debug format in the console and exit \n	-o <origin bound inf> [origin bound sup]\n		set the origin range; if origin bound sup is not given it is set to the most ancient fossil age\n	-e <end bound inf> <end bound sup>\n		set the end time range\n	-p <speciation rate> <extinction rate> <fossilization rate>\n		set the speciation, extinction and fossilization rates\n	-s <number>\n		set the number of samples\n	-d\n		return the distribution (otherwise the density is returned by default)\n	-u <value>\n		set the step discretizing the time distribution to <value>\n	-s <number>\n		set the number of thread running in parallell\n	-h\n		display help\n\n--------------------------\n\n"


static char **readList(FILE *f);



int main(int argc, char **argv) {	
	char *inputFileNameCladeA, *inputFileNameCladeB, *inputFileNameTree, *inputFileNameFossil = NULL, *outputPrefix = PREFIX, outputResult[STRING_SIZE], option[256];
	FILE *fa, *fb, *fi, *ff, *fo;
	int i, j, nSamp = 100, nBurn = 1000, nGap = 10, maxT = 8;
	double minTimeIntInf = NO_TIME, minTimeIntSup = NO_TIME, maxTimeIntInf = 0., step = 0.1, al = 0.75, prop = 0.25, probSpe = 0.33, probExt = 0.33, order = 0.95;
	TypeModelParam param = {.birth=1.3, .death = 1., .fossil = 1., .sampling = 1.}, windSize = {.birth=0.1, .death = 0.1, .fossil = 0.02, .sampling = 1.}, init = {.birth=0.5, .death = 0.5, .fossil = 0.1, .sampling = 1.};

//feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);	
	for(i=0; i<256; i++)
		option[i] = 0;
	   
	for(i=1; i<argc && *(argv[i]) == '-'; i++) {
		for(j=1; argv[i][j] != '\0'; j++)
			option[argv[i][j]] = 1;
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
			} else {
				error("Error while reading %s.\n", argv[i]);
				exit(1);
			}
			exit(0);
		}
		if(option['o']) {
			option['o'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &minTimeIntInf) == 1)
				i++;
			else
				error("at least one value is expected after -o");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &minTimeIntSup) == 1)
				i++;
			else
				minTimeIntSup = NO_TIME;
		}
		if(option['p']) {
			option['p'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(param.birth)) == 1)
				i++;
			else
				error("3 values are expected after -p");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(param.death)) == 1)
				i++;
			else
				error("3 values are expected after -p");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(param.fossil)) == 1)
				i++;
			else
				error("3 values are expected after -p");

		}
		if(option['a']) {
			option['a'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &probSpe) == 1)
				i++;
			else
				error("2 values are expected after -a");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &probExt) == 1)
				i++;
			else
				error("2 values are expected after -a");
		}
		if(option['w']) {
			option['w'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &windSize.birth) == 1)
				i++;
			else
				error("3 values are expected after -w");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &windSize.death) == 1)
				i++;
			else
				error("3 values are expected after -w");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &windSize.fossil) == 1)
				i++;
			else
				error("3 values are expected after -w");
		}
		if(option['i']) {
			option['i'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &init.birth) == 1)
				i++;
			else
				error("3 values are expected after -i");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &init.death) == 1)
				i++;
			else
				error("3 values are expected after -i");
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &init.fossil) == 1)
				i++;
			else
				error("3 values are expected after -i");
		}
		if(option['r']) {
			unsigned int seed;
			option['r'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%u", &seed) == 1) {
				srand(seed);
				i++;
			} else
				error("1 values are expected after -r");
		}
		if(option['s']) {
			option['s'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nSamp) == 1)
				i++;
			else
				error("a number is expected after -s");
		}
		if(option['b']) {
			option['b'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nBurn) == 1)
				i++;
			else
				error("a number is expected after -b");
		}
		if(option['g']) {
			option['g'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &nGap) == 1)
				i++;
			else
				error("a number is expected after -b");
		}
		if(option['l']) {
			option['l'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &al) == 1)
				i++;
			else
				error("a number is expected after -l");
		}
		if(option['u']) {
			option['u'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%le", &(step)) == 1)
				i++;
			else
				error("a number is expected after -u");
		}
		if(option['t']) {
			option['t'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%d", &maxT) == 1)
				i++;
			else
				error("a number is required after option -t");
		}
		if(option['f']) {
			option['f'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &prop) == 1)
				i++;
			else
				error("a number is required after option -t");
		}
		if(option['r']) {
			option['r'] = 0;
			if((i+1)<argc && sscanf(argv[i+1], "%lf", &(order)) == 1)
				i++;
			else
				error("a number is expected after -m");
		}
		if(option['h']) {
			option['h'] = 0;
			printf("%s\n", HELPMESSAGE);
			exit(0);
		}
	}
	if(i<argc) {
		inputFileNameCladeA = argv[i++];
	} else
		error("Please provide the name of a file containing a phylogenetic tree in Newick format\n");
	if(i<argc) {
		inputFileNameCladeB = argv[i++];
	} else
		error("Please provide the name of a file containing a phylogenetic tree in Newick format\n");
	if(i<argc) {
		inputFileNameTree = argv[i++];
	} else
		error("Please provide the name of a file containing a phylogenetic tree in Newick format\n");
	if(i<argc) {
		inputFileNameFossil = argv[i++];
	} else
		error("Please provide the name of a file containing a fossil list\n");
	if(i<argc)
		outputPrefix = argv[i++];
	if((fa = fopen(inputFileNameCladeA, "r")) && (fb = fopen(inputFileNameCladeB, "r")) && (fi = fopen(inputFileNameTree, "r")) && (ff = fopen(inputFileNameFossil, "r"))) {
		TypeTree **tree;
		TypeFossilIntFeature *fos;
		TypeLexiTree *dict;
		int i, n, *cladeA, *cladeB;
		char **list, nameOutput[STRING_SIZE], nameIndex[STRING_SIZE];
		FILE *fout, *find;
		double logProb;

		tree = readTrees(fi);
        fclose(fi);
        if(tree[0] == NULL) {
			error("Error: no tree\n");
			return 1;
		}
		fixTreeSet(tree);
		fos = getFossilIntFeature(ff, tree[0]->name, tree[0]->size);
		fclose(ff);
		if(getMaxFossilIntTime(fos) > 0.)
			negateFossilInt(fos);
		fixStatus(tree[0], fos);
		for(i=1; tree[i]!=NULL; i++) {
			int n;
			for(n=0; n<tree[i]->size; n++)
				if(fos->status[n] == unknownNodeStatus)
					tree[i]->time[n] = fos->endTimeTable[n].inf;
		}

		double minFossilTime = getMinFossilIntTime(fos);
		if(minTimeIntInf == NO_TIME || minTimeIntInf > minFossilTime) {
			if(minFossilTime<0)
				minTimeIntInf = 1.2*minFossilTime;
			else
				minTimeIntInf = 0.8*minFossilTime;
		}
		for(i=0; tree[i]!=NULL; i++) {
			tree[i]->minTime = minTimeIntInf;
			tree[i]->minTimeInt.inf = minTimeIntInf;
			tree[i]->minTimeInt.sup = minTimeIntSup;
			tree[i]->maxTime = maxTimeIntInf;
		}
		sprintf(nameOutput, "%s.out", outputPrefix);
		sprintf(nameIndex, "%s.ind", outputPrefix);
		
		dict = newLexiTree();
		for(n=0; n<tree[0]->size; n++)
			if(tree[0]->name && tree[0]->name[n]) {
				if(addWordLexi(tree[0]->name[n], n, dict)>=0)
					error("Warning! duplicate identifier '%s'\n", tree[0]->name[n]);
			}

		list = readList(fa);
		fclose(fa);
		for(i=0; list[i]!=NULL; i++)
			;
		cladeA = (int*) malloc((i+1)*sizeof(int));
		for(i=0; list[i]!=NULL; i++) {
			cladeA[i] = findWordLexi(list[i], dict);
			if(cladeA[i] < 0)
				error("Error! identifier '%s' not found\n", list[i]);
			else {
				if(!(tree[0]->node[cladeA[i]].child == NOSUCH && fos->status[cladeA[i]] == extinctNodeStatus))
					error("Error! taxon %s is not extinct\n", list[i]);
			}
		}
		cladeA[i] = END_LIST_INT;
		for(i=0; list[i]!=NULL; i++)
			free((void*) list[i]);
		free((void*) list);
		list = readList(fb);
		fclose(fb);
		for(i=0; list[i]!=NULL; i++)
			;
		cladeB = (int*) malloc((i+1)*sizeof(int));
		for(i=0; list[i]!=NULL; i++) {
			cladeB[i] = findWordLexi(list[i], dict);
			if(cladeB[i] < 0)
				error("Error! identifier '%s' not found\n", list[i]);
			else {
				if(!(tree[0]->node[cladeB[i]].child == NOSUCH && fos->status[cladeB[i]] == extinctNodeStatus))
					error("Error! taxon %s is not extinct\n", list[i]);
			}
		}
		cladeB[i] = END_LIST_INT;
		freeLexiTree(dict);
		for(i=0; list[i]!=NULL; i++)
			free((void*) list[i]);
		free((void*) list);
		if((fout = fopen(nameOutput, "w")) && (find = fopen(nameIndex, "w"))) {
			logProb = MCMCSamplingExtinctionComp(fout, find, cladeA, cladeB, tree, fos, al, nBurn, nGap, nSamp, prop, &windSize, &init, probSpe, probExt);
			fclose(fout);
			fclose(find);
		} else
			error("Error. Can't write coda files\n");
		sprintf(outputResult, "%s.txt", outputPrefix);
		if((fo = fopen(outputResult, "w"))) {
			fprintf(fo, "probability %lf\n\n", exp(logProb));
			fclose(fo);
		}
fprintf(stdout, "probability %lf\n\n", exp(logProb));
		freeFossilIntFeature(fos);
		for(i=0; tree[i]!=NULL; i++)
			freeTree(tree[i]);
		free((void*)tree);
		free((void*)cladeA);
		free((void*)cladeB);
	} else {
		error("Cannot read %s or %s or %s or %s\n", inputFileNameCladeA, inputFileNameCladeB, inputFileNameTree, inputFileNameFossil);
		exit(1);
	}
	return 0;
}


#define MAX_SIZE_TMP 50
#define INC_BUFFER 50
#define IS_SEP(c) (c == ' ' || c == '\t' || c == '\r' || c == '\n' || c == ';')

char **readList(FILE *f) {
	char c, tmp[MAX_SIZE_TMP+1], **list;
	int size, sizeBuffer;

	sizeBuffer = INC_BUFFER;
	list= (char**) malloc(sizeBuffer*sizeof(char*));
	size = 0;
	do {
		c = getc(f);
	} while(c!=EOF && IS_SEP(c)); 
	while(c != EOF) {
		int i;
		i = 0;
		while(i<MAX_SIZE_TMP && c !=EOF && !IS_SEP(c)) {
			tmp[i] = c;
			c = getc(f);
			i++;
		}
		tmp[i++] = '\0';
		if(i == MAX_SIZE_TMP)
			error("Ident too long (%s) ...", tmp);
		if(i>1) {
			if(size >= sizeBuffer) {
				sizeBuffer += INC_BUFFER;
				list = (char**) realloc((void *) list, sizeBuffer*sizeof(char*));
			}
			list[size] = (char*) malloc((strlen(tmp)+1)*sizeof(char));
			strcpy(list[size], tmp);
			size++;
		}
		while(c!=EOF && IS_SEP(c))
			c=getc(f);
	}
	if(size >= sizeBuffer) {
		sizeBuffer += INC_BUFFER;
		list = (char**) realloc((void *) list, sizeBuffer*sizeof(char*));
	}
	list[size++] = NULL;
	return list;
}
