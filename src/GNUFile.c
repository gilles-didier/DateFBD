#include "GNUFile.h"


void fprintGNUFile(FILE *f, TypeGNUInfo *info) {
	int c;
	if(info->dataFileName == NULL || info->nColumn == 0)
		return;
	switch(info->type) {
		case EPS:
		default:
			fprintf(f, "set terminal postscript enhanced color eps\n");
	}
	if(info->xLabel != NULL)
		fprintf(f, "set xlabel \"%s\"\n",  info->xLabel);
	if(info->yLabel != NULL)
		fprintf(f, "set ylabel \"%s\"\n",  info->yLabel);
	if(info->outputFileName != NULL)
		fprintf(f, "set output \"%s\"\n", info->outputFileName);
	fprintf(f, "plot \"%s\"  using 1:2 with lines", info->dataFileName);
	if(info->title != NULL && info->title[0] != NULL)
		fprintf(f, " title '%s'", info->title[0]);
	for(c=1; c<info->nColumn; c++) {
		fprintf(f,", '' using 1:%d with lines", c+2);
		if(info->title != NULL && info->title[c] != NULL)
			fprintf(f," title '%s'", info->title[c]);
	}
}
