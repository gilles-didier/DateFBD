#ifndef GNUFileF
#define GNUFileF

#include <stdio.h>

typedef enum TYPE_GNU_FILE {
	EPS=0
} TypeGNUFile;


typedef struct GNUFILE_INFO {
	char *outputFileName, *xLabel, *yLabel, *dataFileName, **title;
	int nColumn;
	TypeGNUFile type;
} TypeGNUInfo;

void fprintGNUFile(FILE *f, TypeGNUInfo *info);
	
#endif
