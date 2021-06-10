#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "Utils.h"
#include "FossilInt.h"

#define INC_FOSSIL_ITEM 50
#define MAX_NAME_SIZE 3000

/*add a new fossil to n*/
void addFossilInt(TypeFossilIntFeature *fos, int n, double inf, double sup) {
    if(fos == NULL)
        return;
    if(fos->sizeFossil>=fos->sizeBufFossil) {
        fos->sizeBufFossil += INC_FOSSIL_ITEM;
        fos->fossilIntList = (TypeFossilIntList*) realloc((void*) fos->fossilIntList, fos->sizeBufFossil*sizeof(TypeFossilIntList));
    }
    fos->fossilIntList[fos->sizeFossil].fossilInt.inf = inf;
    fos->fossilIntList[fos->sizeFossil].fossilInt.sup = sup;
    fos->fossilIntList[fos->sizeFossil].prec = fos->fossilInt[n];
    fos->fossilInt[n] = fos->sizeFossil;
    fos->sizeFossil++;
}

/*allocate a new fossil int feature with n entries/nodes*/
TypeFossilIntFeature *newFossilIntFeature(int n) {
    TypeFossilIntFeature *res;
    int i;
    if(n == 0)
        return NULL;
    res = (TypeFossilIntFeature*) malloc(sizeof(TypeFossilIntFeature));
    res->sizeFossil = 0;
    res->sizeBufFossil = INC_FOSSIL_ITEM;
    res->sizeNode = 0;
    res->sizeBufNode = n;
    res->sizeTime = 0;
    res->sizeBufTime = n/2;
    res->fossilInt = (int*) malloc(res->sizeBufNode*sizeof(int));
    for(i=0; i<res->sizeBufNode; i++)
        res->fossilInt[i] = NOSUCH;
    res->fossilIntList = (TypeFossilIntList*) malloc(res->sizeBufFossil*sizeof(TypeFossilIntList));
    res->status = (TypeNodeStatus*) malloc(res->sizeBufNode*sizeof(TypeNodeStatus));
    for(i=0; i<res->sizeBufNode; i++)
        res->status[i] = noneNodeStatus;
    res->endTime = (int*) malloc(res->sizeBufNode*sizeof(int));
    for(i=0; i<res->sizeBufNode; i++)
        res->endTime[i] = NOSUCH;
    res->endTimeTable = (TypeTimeInterval*) malloc(res->sizeBufTime*sizeof(TypeTimeInterval));
    return res;
}

/*fix status of leaves*/
void fixStatus(TypeTree *tree, TypeFossilIntFeature *feat) {
    int n;
    for(n=0; n<tree->size; n++) {
		if(feat->status[n] == unknownNodeStatus)
			tree->time[n] = feat->endTimeTable[n].inf;
        if(tree->node[n].child == NOSUCH && feat->status[n] == noneNodeStatus) {
            if(tree->time == NULL || tree->time[n] == NO_TIME) {
                if(feat->fossilInt[n] == NOSUCH)
                    feat->status[n] = contempNodeStatus;
                else
                    feat->status[n] = extinctNodeStatus;
            } else
                feat->status[n] = unknownNodeStatus;
        }
    }
}

/*fully duplicate "feat"*/
TypeFossilIntFeature *cpyFossilIntFeature(TypeFossilIntFeature *feat, int n) {
    TypeFossilIntFeature *res;
    int i;
    if(feat == NULL)
        return NULL;
    res = (TypeFossilIntFeature*) malloc(sizeof(TypeFossilIntFeature));
    res->sizeFossil = feat->sizeFossil;
    res->sizeBufFossil = feat->sizeFossil;
    res->fossilInt = (int*) malloc(n*sizeof(int));
    for(i=0; i<n; i++)
        res->fossilInt[i] = feat->fossilInt[i];
    res->fossilIntList = (TypeFossilIntList*) malloc(res->sizeBufFossil*sizeof(TypeFossilIntList));
    for(i=0; i<res->sizeFossil; i++)
        res->fossilIntList[i] = feat->fossilIntList[i];
    res->sizeBufNode = feat->sizeBufNode;
    res->status = (TypeNodeStatus*) malloc(res->sizeBufNode*sizeof(TypeNodeStatus));
    for(i=0; i<res->sizeBufNode; i++)
        res->status[i] = feat->status[i];
    res->endTime = (int*) malloc(res->sizeBufNode*sizeof(int));
    for(i=0; i<res->sizeBufNode; i++)
        res->endTime[i] = feat->endTime[i];
    res->endTimeTable = (TypeTimeInterval*) malloc(res->sizeBufTime*sizeof(TypeTimeInterval));
    for(i=0; i<res->sizeTime; i++)
        res->endTimeTable[i] = feat->endTimeTable[i];
    return res;
}


void fixTree(TypeTree *tree, TypeFossilIntFeature *fos) {
    int n;
    if(fos == NULL || tree->time == NULL)
        return;
    for(n=0; n<tree->size; n++)
        if(tree->node[n].child == NOSUCH && tree->time[n] == NO_TIME && fos->fossilInt[n] != NOSUCH) {
            int f;
            tree->time[n] = fos->fossilIntList[fos->fossilInt[n]].fossilInt.sup;
            for(f=fos->fossilIntList[fos->fossilInt[n]].prec; f!=NOSUCH; f=fos->fossilIntList[f].prec)
                if(fos->fossilIntList[f].fossilInt.sup>tree->time[n])
                    tree->time[n] = fos->fossilIntList[f].fossilInt.sup;
        }
}

double getMinFossilIntTimeFromNode(int n, TypeTree *tree, TypeFossilIntFeature* feat) {
	double min = tree->maxTime;
	int f;
	if(feat == NULL || feat->sizeFossil == 0)
		return 0.;
	if(tree->time[n] != NO_TIME)
		min = tree->time[n];
	if(feat->fossilInt[n] != NOSUCH) {
		min = feat->fossilIntList[feat->fossilInt[n]].fossilInt.sup;
		for(f=feat->fossilIntList[feat->fossilInt[n]].prec; f!=NOSUCH; f=feat->fossilIntList[f].prec)
			if(feat->fossilIntList[f].fossilInt.sup<min)
				min = feat->fossilIntList[f].fossilInt.sup;
	} else {
		if(tree->node[n].child != NOSUCH) {
			int c;
			min = getMinFossilIntTimeFromNode(tree->node[n].child, tree, feat);
			for(c=tree->node[tree->node[n].child].sibling; c!=NOSUCH; c=tree->node[c].sibling) {
				double tmp = getMinFossilIntTimeFromNode(c, tree, feat);
				if(tmp<min)
					min = tmp;
			}
		}
	}
	return min;
}

double getMaxFossilIntTimeToNode(int n, TypeTree *tree, TypeFossilIntFeature* feat) {
	double max = tree->minTimeInt.sup;
	if(feat->fossilInt[n] != NOSUCH) {
		int f;
		max = feat->fossilIntList[feat->fossilInt[n]].fossilInt.inf;
		for(f=feat->fossilIntList[feat->fossilInt[n]].prec; f!=NOSUCH; f=feat->fossilIntList[f].prec)
			if(feat->fossilIntList[f].fossilInt.inf>max)
				max = feat->fossilIntList[f].fossilInt.inf;
	} else {
		if(tree->parent[n] != NOSUCH)
			max = getMaxFossilIntTimeToNode(tree->parent[n], tree, feat);
		else
			max = tree->minTimeInt.inf;
	}
	return max;
}

double getFirstFossilIntTimeFromNode(int n, TypeTree *tree, TypeFossilIntFeature* feat) {
	double min = DBL_MAX;
	if(feat == NULL || feat->sizeFossil == 0)
		return 0.;
	if(feat->fossilInt[n] != NOSUCH) {
		int f;
		min = feat->fossilIntList[feat->fossilInt[n]].fossilInt.sup;
		for(f=feat->fossilIntList[feat->fossilInt[n]].prec; f!=NOSUCH; f=feat->fossilIntList[f].prec)
			if(feat->fossilIntList[f].fossilInt.sup<min)
				min = feat->fossilIntList[f].fossilInt.sup;
	} else {
		if(tree->node[n].child != NOSUCH) {
			int c;
			min = getFirstFossilIntTimeFromNode(tree->node[n].child, tree, feat);
			for(c=tree->node[tree->node[n].child].sibling; c!=NOSUCH; c=tree->node[c].sibling) {
				double tmp = getFirstFossilIntTimeFromNode(c, tree, feat);
				if(tmp<min)
				min = tmp;
			}
		}
	}
	return min;
}
	
double getMidFossilIntTimeFromNode(int n, TypeTree *tree, TypeFossilIntFeature* feat) {
	double min = DBL_MAX;
	if(feat == NULL || feat->sizeFossil == 0)
		return 0.;
	if(feat->fossilInt[n] != NOSUCH) {
		int f;
		min = (feat->fossilIntList[feat->fossilInt[n]].fossilInt.sup+feat->fossilIntList[feat->fossilInt[n]].fossilInt.inf)/2.;
		for(f=feat->fossilIntList[feat->fossilInt[n]].prec; f!=NOSUCH; f=feat->fossilIntList[f].prec) {
			double tmp = (feat->fossilIntList[f].fossilInt.sup+feat->fossilIntList[f].fossilInt.inf)/2.;
			if(tmp<min)
				min = tmp;	
		}
	} else {
		if(tree->node[n].child != NOSUCH) {
			int c;
			min = getMidFossilIntTimeFromNode(tree->node[n].child, tree, feat);
			for(c=tree->node[tree->node[n].child].sibling; c!=NOSUCH; c=tree->node[c].sibling) {
				double tmp = getMidFossilIntTimeFromNode(c, tree, feat);
				if(tmp<min)
				min = tmp;
			}
		}
	}
	return min;
}


double getFirstFossilIntTime(TypeFossilIntFeature* feat) {
    int i;
    double min;
    if(feat == NULL || feat->sizeFossil == 0)
        return 0.;
    min = feat->fossilIntList[0].fossilInt.sup;
    for(i=1; i<feat->sizeFossil; i++)
        if(feat->fossilIntList[i].fossilInt.sup<min)
            min = feat->fossilIntList[i].fossilInt.sup;
    return min;
}

double getMidFossilIntTime(TypeFossilIntFeature* feat) {
    int i;
    double min;
    if(feat == NULL || feat->sizeFossil == 0)
        return 0.;
    min = (feat->fossilIntList[0].fossilInt.inf+feat->fossilIntList[0].fossilInt.sup)/2.;
    for(i=1; i<feat->sizeFossil; i++)
        if((feat->fossilIntList[i].fossilInt.inf+feat->fossilIntList[i].fossilInt.sup)/2.<min)
            min = (feat->fossilIntList[i].fossilInt.inf+feat->fossilIntList[i].fossilInt.sup)/2.;
    return min;
}

double getMinFossilIntTime(TypeFossilIntFeature* feat) {
    int i;
    double min;
    if(feat == NULL || feat->sizeFossil == 0)
        return 0.;
    min = feat->fossilIntList[0].fossilInt.inf;
    for(i=1; i<feat->sizeFossil; i++)
        if(feat->fossilIntList[i].fossilInt.inf<min)
            min = feat->fossilIntList[i].fossilInt.inf;
    return min;
}

double getMaxFossilIntTime(TypeFossilIntFeature* feat) {
    int i;
    double max;
    if(feat == NULL || feat->sizeFossil == 0)
        return 0.;
    max = feat->fossilIntList[0].fossilInt.sup;
    for(i=1; i<feat->sizeFossil; i++)
        if(feat->fossilIntList[i].fossilInt.sup>max)
            max = feat->fossilIntList[i].fossilInt.sup;
    return max;
}

double getFossilIntMaxInfTime(TypeFossilIntFeature* feat) {
    int i;
    double max;
    if(feat == NULL || feat->sizeFossil == 0)
        return 0.;
    max = feat->fossilIntList[0].fossilInt.inf;
    for(i=1; i<feat->sizeFossil; i++)
        if(feat->fossilIntList[i].fossilInt.inf>max)
            max = feat->fossilIntList[i].fossilInt.inf;
    return max;
}

void negateFossilInt(TypeFossilIntFeature* feat) {
    int i;
    for(i=0; i<feat->sizeFossil; i++) {
        double tmp = feat->fossilIntList[i].fossilInt.inf;
        feat->fossilIntList[i].fossilInt.inf = -feat->fossilIntList[i].fossilInt.sup;
        feat->fossilIntList[i].fossilInt.sup = -tmp;
    }
    if(feat->endTimeTable != NULL) {
		for(i=0; i<feat->sizeNode; i++) {
			double tmp = feat->endTimeTable[i].inf;
			feat->endTimeTable[i].inf = -feat->endTimeTable[i].sup;
			feat->endTimeTable[i].sup = -tmp;
		}
    }		
}



TypeFossilFeature *sampleFossilInt(TypeFossilIntFeature* feat, int size) {
    TypeFossilFeature *sample;
    int i, n, f, *tmp;
    size_t *index;

    if(feat == NULL)
        return NULL;
    sample = (TypeFossilFeature*) malloc(sizeof(TypeFossilFeature));
    sample->size = feat->sizeFossil;
    sample->sizeBuf = feat->sizeFossil;
    sample->fossilList = (TypeFossilList*) malloc(sample->size*sizeof(TypeFossilList));
    sample->fossil = (int*) malloc(size*sizeof(int));
    sample->status = (TypeNodeStatus*) malloc(size*sizeof(int));
    for(i=0; i<size; i++)
        sample->fossil[i] = feat->fossilInt[i];
    for(i=0; i<sample->size; i++) {
        sample->fossilList[i].prec = feat->fossilIntList[i].prec;
        sample->fossilList[i].time = feat->fossilIntList[i].fossilInt.inf+UNIF_RAND*(feat->fossilIntList[i].fossilInt.sup-feat->fossilIntList[i].fossilInt.inf);
    }
    index = qsortindex(sample->fossilList, sample->size, sizeof(TypeFossilList), compareFossilList);
    for(n=0; n<size; n++)
        if(sample->fossil[n]>=0)
            sample->fossil[n] = index[sample->fossil[n]];
    for(f=0; f<sample->size; f++)
        if(sample->fossilList[f].prec>=0)
            sample->fossilList[f].prec = index[sample->fossilList[f].prec];
    free((void*) index);
    tmp = (int*) malloc((sample->size+1)*sizeof(int));
    for(n=0; n<size; n++) {
        int ind = 0;
        for(f=sample->fossil[n]; f>=0; f=sample->fossilList[f].prec)
            tmp[ind++] = f;
        if(ind>0) {
            qsort(tmp, ind, sizeof(int), compareInt);
            sample->fossilList[tmp[0]].prec = -1;
            for(f=1; f<ind; f++)
                sample->fossilList[tmp[f]].prec = tmp[f-1];
            sample->fossil[n] = tmp[ind-1];
        }
    }
    free((void*) tmp);
    return sample;
}


TypeFossilFeature *setTimeFossilInt(TypeFossilIntFeature* feat, double *time, int size) {
    TypeFossilFeature *sample;
    int i, n, f, *tmp;
    size_t *index;

    if(feat == NULL)
        return NULL;
    sample = (TypeFossilFeature*) malloc(sizeof(TypeFossilFeature));
    sample->size = feat->sizeFossil;
    sample->sizeBuf = feat->sizeFossil;
    sample->fossilList = (TypeFossilList*) malloc(sample->size*sizeof(TypeFossilList));
    sample->fossil = (int*) malloc(size*sizeof(int));
    sample->status = (TypeNodeStatus*) malloc(size*sizeof(int));
    for(i=0; i<size; i++)
        sample->fossil[i] = feat->fossilInt[i];
    for(i=0; i<sample->size; i++) {
        sample->fossilList[i].prec = feat->fossilIntList[i].prec;
		sample->fossilList[i].time = time[i];
    }
    index = qsortindex(sample->fossilList, sample->size, sizeof(TypeFossilList), compareFossilList);
    for(n=0; n<size; n++)
        if(sample->fossil[n]>=0)
            sample->fossil[n] = index[sample->fossil[n]];
    for(f=0; f<sample->size; f++)
        if(sample->fossilList[f].prec>=0)
            sample->fossilList[f].prec = index[sample->fossilList[f].prec];
    free((void*) index);
    tmp = (int*) malloc((sample->size+1)*sizeof(int));
    for(n=0; n<size; n++) {
        int ind = 0;
        for(f=sample->fossil[n]; f>=0; f=sample->fossilList[f].prec)
            tmp[ind++] = f;
        if(ind>0) {
            qsort(tmp, ind, sizeof(int), compareInt);
            sample->fossilList[tmp[0]].prec = -1;
            for(f=1; f<ind; f++)
                sample->fossilList[tmp[f]].prec = tmp[f-1];
            sample->fossil[n] = tmp[ind-1];
        }
    }
    free((void*) tmp);
    return sample;
}


/*fully duplicate "tab"*/
TypeFossilIntTab cpyFossilIntTab(TypeFossilIntTab tab) {
    TypeFossilIntTab res;
    int i;
    res.size = tab.size;
    res.fossilInt = (TypeTimeInterval*) malloc(res.size*sizeof(TypeTimeInterval));
    for(i=0; i<res.size; i++)
        res.fossilInt[i] = tab.fossilInt[i];
    return res;
}


TypeTimeInterval toFossilInt(char *s) {
    int i, ind;
    TypeTimeInterval fossilInt;
    char tmp[MAX_NAME_SIZE];
    for(i=0; s[i] != '\0'; i++)
		if(s[i] == ',')
			s[i] = '.';
    for(i=0; s[i] != '\0' && isspace(s[i]); i++);
    if(s[i] == '(' || s[i] == '[') {
        i++;
        for(; s[i] != '\0' && isspace(s[i]); i++);
        ind = 0;
        for(; s[i] != '\0' && s[i] != ':' && s[i] != ';' && !isspace(s[i]); i++)
            tmp[ind++] = s[i];
        tmp[ind] = '\0';
        fossilInt.inf =	atof(tmp);
        for(; s[i] != '\0' && isspace(s[i]); i++);
        if(s[i] != ':' && s[i] != ';')
            error("Missing ':' or ';' while reading a fossilInt time interval.");
        i++;
        for(; s[i] != '\0' && isspace(s[i]); i++);
        ind = 0;
        for(; s[i] != '\0' && s[i] != ')' && s[i] != ']' && !isspace(s[i]); i++)
            tmp[ind++] = s[i];
        tmp[ind] = '\0';
        fossilInt.sup =	atof(tmp);
        for(; s[i] != '\0' && isspace(s[i]); i++);
        if(s[i] != ')' && s[i] != ']')
            error("Missing ')' or ']' while reading a fossilInt time interval.");
        if(fossilInt.sup<fossilInt.inf) {
            double tmp = fossilInt.sup;
            fossilInt.sup = fossilInt.inf;
            fossilInt.inf = tmp;
        }
    } else {
        fossilInt.inf = atof(s);
        fossilInt.sup = fossilInt.inf;
    }
    return fossilInt;
}

int compareFossilInt(const void* a, const void* b) {
    if(((TypeTimeInterval*)a)->inf>((TypeTimeInterval*)b)->inf)
        return 1;
    if(((TypeTimeInterval*)a)->inf<((TypeTimeInterval*)b)->inf)
        return -1;
    if(((TypeTimeInterval*)a)->sup>((TypeTimeInterval*)b)->sup)
        return 1;
    if(((TypeTimeInterval*)a)->sup<((TypeTimeInterval*)b)->sup)
        return -1;
    return 0;
}

int compareFossilIntList(const void* a, const void* b) {
    if(((TypeFossilIntList*)a)->fossilInt.inf>((TypeFossilIntList*)b)->fossilInt.inf)
        return 1;
    if(((TypeFossilIntList*)a)->fossilInt.inf<((TypeFossilIntList*)b)->fossilInt.inf)
        return -1;
    if(((TypeFossilIntList*)a)->fossilInt.sup>((TypeFossilIntList*)b)->fossilInt.sup)
        return 1;
    if(((TypeFossilIntList*)a)->fossilInt.sup<((TypeFossilIntList*)b)->fossilInt.sup)
        return -1;
    return 0;
}

/*get the fossilInts from comments*/
TypeFossilIntTab getFossilInt(char *str) {
    int i, ind = 0, sizeBufFossil = INC_FOSSIL_ITEM;
    char tmp[MAX_NAME_SIZE+1];
    TypeFossilIntTab tab;
    tab.size = 0;
    if(str == NULL || strlen(str) < 5) {
        tab.fossilInt = NULL;
        return tab;
    } else
        tab.fossilInt = (TypeTimeInterval*) malloc(sizeBufFossil*sizeof(TypeTimeInterval));
    i = 4;
    while(str[i] != '\0') {
        for(; str[i] != '\0' && (str[i-4]!=':' || str[i-3]!='F' || str[i-2]!='O' || str[i-1]!='S' || str[i]!='='); i++);
        if(str[i]=='=')
            i++;
        if(str[i] != '\0') {
            for(; str[i] != '\0' && issep(str[i]); i++);
            ind = 0;
            for(; str[i] != '\0' && !issep(str[i]) && ind<MAX_NAME_SIZE; i++)
                tmp[ind++] = str[i];
            tmp[ind++] = '\0';
            if(tab.size>=sizeBufFossil) {
                sizeBufFossil += INC_FOSSIL_ITEM;
                tab.fossilInt = (TypeTimeInterval*) realloc((void*) tab.fossilInt, sizeBufFossil*sizeof(TypeTimeInterval));
            }
            tab.fossilInt[tab.size++] = toFossilInt(tmp);
        }
    }
    if(tab.size>0) {
        tab.fossilInt = (TypeTimeInterval*) realloc((void*) tab.fossilInt, tab.size*sizeof(TypeTimeInterval));
        for(i=0; i<tab.size; i++)
            if(tab.fossilInt[i].inf>tab.fossilInt[i].sup) {
                double tmp = tab.fossilInt[i].inf;
                tab.fossilInt[i].inf = tab.fossilInt[i].sup;
                tab.fossilInt[i].sup = tmp;
            }
        qsort(tab.fossilInt, tab.size, sizeof(TypeTimeInterval), compareFossilInt);
    } else {
        free((void*) tab.fossilInt);
        tab.fossilInt = NULL;
    }
    return tab;
}


/*from fos to FossilInt fos*/
TypeFossilIntFeature *fosToFossilInt(TypeTree *tree) {
    TypeFossilIntFeature *fos;
    int n;
    char *s;
    fos = newFossilIntFeature(tree->sizeBuf);
    fos->sizeNode = tree->size;
    if(tree->comment != NULL) {
        for(n=0; n<tree->size; n++) {
           if(tree->comment[n] != NULL) {
                TypeTimeInterval ti;
                ti.inf = sqrt(-1);
                ti.sup = sqrt(-1);
                s = tree->comment[n];
                while(s!=NULL && s[0]!='\0') {
                    char *tag, *val;
                    s = nextTag(s, &tag, &val);
                    if(tag != NULL && strcmp(tag, "FOS")==0 && val != NULL) {
                        if(fos->sizeFossil>=fos->sizeBufFossil) {
                            fos->sizeBufFossil += INC_FOSSIL_ITEM;
                            fos->fossilIntList = (TypeFossilIntList*) realloc((void*) fos->fossilIntList, fos->sizeBufFossil*sizeof(TypeFossilIntList));
                        }
                        fos->fossilIntList[fos->sizeFossil].fossilInt = toFossilInt(val);
                        fos->fossilIntList[fos->sizeFossil].prec = fos->fossilInt[n];
                        fos->fossilInt[n] = fos->sizeFossil++;
                    }
                    if(tag != NULL && strcmp(tag, TAG_STATUS)==0 && val != NULL) {
                        int st = atoi(val);
                        switch(st) {
                            case 1:
                                fos->status[n] = contempNodeStatus;
                            break;
                            case 2:
                                fos->status[n] = extinctNodeStatus;
                            break;
                            case 3:
                                fos->status[n] = unknownNodeStatus;
                            break;
                            default:
                                fos->status[n] = noneNodeStatus;
                        }
                    }
                    if(tag != NULL && strcmp(tag, TAG_TIME)==0 && val != NULL)
                        ti = toFossilInt(val);
                }
               if(fos->status[n] == unknownNodeStatus) {
                    if(!isnan(ti.inf)) {
                         tree->time[n] = ti.inf;
                   } else {
                        if(tree->node[n].child == NOSUCH)
                            fos->status[n] = contempNodeStatus;
                        else
                            fos->status[n] = noneNodeStatus;
                    }
                }
            }
        }
        if(tree->comment[tree->root] != NULL) {
            TypeTimeInterval ti;
            ti.inf = sqrt(-1);
            ti.sup = sqrt(-1);
            s = tree->comment[tree->root];
            while(s!=NULL && s[0]!='\0') {
                char *tag, *val;
                s = nextTag(s, &tag, &val);
                if(tag != NULL && strcmp(tag, TAG_FOSSIL)==0 && val != NULL) {
                    if(fos->sizeFossil>=fos->sizeBufFossil) {
                        fos->sizeBufFossil += INC_FOSSIL_ITEM;
                        fos->fossilIntList = (TypeFossilIntList*) realloc((void*) fos->fossilIntList, fos->sizeBufFossil*sizeof(TypeFossilIntList));
                    }
                    fos->fossilIntList[fos->sizeFossil].fossilInt = toFossilInt(val);
                    fos->fossilIntList[fos->sizeFossil].prec = fos->fossilInt[tree->root];
                    fos->fossilInt[tree->root] = fos->sizeFossil++;
                }
                if(tag != NULL && strcmp(tag, TAG_ORIGIN)==0 && val != NULL)
                    tree->minTimeInt = toFossilInt(val);
                if(tag != NULL && strcmp(tag, TAG_END)==0 && val != NULL)
                    tree->maxTimeInt = toFossilInt(val);
                if(tag != NULL && strcmp(tag, TAG_STATUS)==0 && val != NULL) {
                    int st = atoi(val);
                    switch(st) {
                        case 1:
                            fos->status[tree->root] = contempNodeStatus;
                        break;
                        case 2:
                            fos->status[tree->root] = extinctNodeStatus;
                        break;
                        case 3:
                            fos->status[tree->root] = unknownNodeStatus;
                        break;
                        default:
                            fos->status[tree->root] = noneNodeStatus;
                    }
                }
                if(tag != NULL && strcmp(tag, TAG_TIME)==0 && val != NULL)
                    ti = toFossilInt(val);
            }
            if(fos->status[tree->root] == unknownNodeStatus) {
                if(!isnan(ti.inf)) {
                    if(fos->sizeTime>=fos->sizeBufTime) {
                        fos->sizeBufTime += INC_BUF_TIME_TABLE;
                        fos->endTimeTable = (TypeTimeInterval*) realloc ((void*)fos->endTimeTable,fos->sizeBufTime*sizeof(TypeTimeInterval));
                    }
                    fos->endTimeTable[fos->sizeTime] = ti;
                    fos->endTime[tree->root] = fos->sizeTime++;
                } else {
                    if(tree->node[tree->root].child == NOSUCH)
                        fos->status[tree->root] = contempNodeStatus;
                    else
                        fos->status[tree->root] = noneNodeStatus;
                }
            }
        }
    }
    fixStatus(tree, fos);
    return fos;
}

TypeNameFossilIntTab *readFossilIntTab(FILE *f) {
    int sizeBuf;
    char c;
    TypeNameFossilIntTab *res;

    res = (TypeNameFossilIntTab *) malloc(sizeof(TypeNameFossilIntTab));
    sizeBuf = INC_FOSSIL_ITEM;
    res->name = (char**) malloc(sizeBuf*sizeof(char*));
    res->stopTime = (double*) malloc(sizeBuf*sizeof(double));
    res->fossilIntTab = (TypeFossilIntTab*) malloc(sizeBuf*sizeof(TypeFossilIntTab));
    res->size = 0;
    do {
        char *tmp;
        int i;
        tmp = (char*) malloc((MAX_NAME_SIZE+1)*sizeof(char));
        for(c = fgetc(f); c != EOF && issepline(c); c = fgetc(f));
        if(c == '\'' || c == '"') {
            c = fgetc(f);
            for(i=0; i<MAX_NAME_SIZE && c != EOF && c != '\'' && c != '"'; i++) {
                tmp[i] = c;
                c = fgetc(f);
            }
            if(c == '\'' || c == '"')
                c = fgetc(f);
        } else {
            for(i=0; i<MAX_NAME_SIZE && c !=EOF && !issep(c); i++) {
                tmp[i] = c;
                c = fgetc(f);
            }
        }
        if(i == MAX_NAME_SIZE) {
			tmp[i] = '\0';
            error("Name too much long:\n'%s'", tmp);
        }
        tmp[i++] = '\0';
        if(i>1) {
            char bof[MAX_NAME_SIZE+1];
            int sizeFossilIntBuf = INC_FOSSIL_ITEM;
            if(res->size >= sizeBuf) {
                sizeBuf += INC_FOSSIL_ITEM;
                res->name = (char**) realloc((void*) res->name, sizeBuf*sizeof(char*));
                res->stopTime = (double*) realloc((void*) res->stopTime, sizeBuf*sizeof(double));
				res->fossilIntTab = (TypeFossilIntTab*) realloc((void*) res->fossilIntTab, sizeBuf*sizeof(TypeFossilIntTab));
            }
            res->name[res->size] = (char *) realloc((void*) tmp, i*sizeof(char));
            res->stopTime[res->size] = NO_TIME;
            fixSpace(res->name[res->size]);
            res->fossilIntTab[res->size].size = 0;
            res->fossilIntTab[res->size].fossilInt = (TypeTimeInterval*) malloc(sizeFossilIntBuf*sizeof(TypeTimeInterval));
            for(; c != EOF && issep(c); c = fgetc(f));
            while(c != EOF && !isline(c)) {
                for(i=0; c != EOF && !issepline(c) && i<MAX_NAME_SIZE; i++) {
                    bof[i] = c;
                    c = fgetc(f);
                }
                bof[i++] = '\0';
                if(res->fossilIntTab[res->size].size>=sizeFossilIntBuf) {
                    sizeFossilIntBuf += INC_FOSSIL_ITEM;
                    res->fossilIntTab[res->size].fossilInt = (TypeTimeInterval*) realloc((void*) res->fossilIntTab[res->size].fossilInt, sizeFossilIntBuf*sizeof(TypeTimeInterval));
                }
                if(bof[0] == 'S' || bof[0] == 's') {
					res->stopTime[res->size] = atof(&(bof[1]));
				} else {
					res->fossilIntTab[res->size].fossilInt[res->fossilIntTab[res->size].size++] = toFossilInt(bof);
				}
               for(; c != EOF && issep(c); c = fgetc(f));
            }
            res->size++;
        } else
            free((void*) tmp);
    } while(c != EOF);
    if(res->size > 0)
        res->fossilIntTab = (TypeFossilIntTab*) realloc((void*) res->fossilIntTab, (res->size)*sizeof(TypeFossilIntTab));
    else {
        free((void*)res->fossilIntTab);
        res->fossilIntTab = NULL;
    }
    return res;
}

void freeNameFossilIntTab(TypeNameFossilIntTab *list) {
    int i;
    for(i=0; i<list->size; i++) {
        if(list->name[i] != NULL)
            free((void*)list->name[i]);
        if(list->fossilIntTab[i].fossilInt != NULL)
            free((void*)list->fossilIntTab[i].fossilInt);
    }
    free((void*)list->name);
    free((void*)list->fossilIntTab);
    free((void*)list);
}

void freeIndexFossilIntTab(TypeIndexFossilIntTab *list) {
    int i;
    for(i=0; i<list->size; i++)
        if(list->fossilIntTab[i].fossilInt != NULL)
            free((void*)list->fossilIntTab[i].fossilInt);
    free((void*)list->index);
    free((void*)list->fossilIntTab);
    free((void*)list);
}

   TypeFossilIntList *fossilIntList;
    TypeTimeInterval *endTimeTable;
    int *fossilInt, sizeFossil, sizeBufFossil, sizeTime, sizeBufTime, sizeNode, sizeBufNode, *endTime;
    TypeNodeStatus *status;


void freeFossilIntFeature(TypeFossilIntFeature *fos) {
	if(fos == NULL)
		return;
	if(fos->fossilIntList != NULL)
		free((void*)fos->fossilIntList);
	if(fos->endTimeTable != NULL)
		free((void*)fos->endTimeTable);
	if(fos->fossilInt != NULL)
		free((void*)fos->fossilInt);
	if(fos->endTime != NULL)
		free((void*)fos->endTime);
	if(fos->status != NULL)
		free((void*)fos->status);
	free((void*)fos);
}

void removeFossilIntFossil(TypeFossilIntFeature *fos, int* list) {
    if(fos == NULL || list == NULL)
        return;
    int size = 0, f, n, *index = (int*) malloc(fos->sizeFossil*sizeof(int));
    for(f=0; f<fos->sizeFossil; f++)
        index[f] = 0;
    for(f=0; list[f]!=NOSUCH; f++)
        index[list[f]] = NOSUCH;
    for(f=0; f<fos->sizeFossil; f++)
        if(index[f] != NOSUCH)
            index[f] = size++;
    for(n=0; n<fos->sizeNode; n++)
        if(fos->fossilInt[n] != NOSUCH)
            fos->fossilInt[n] = index[fos->fossilInt[n]];
    for(f=0; f<fos->sizeFossil; f++)
        if(index[f] != NOSUCH) {
            fos->fossilIntList[index[f]].fossilInt = fos->fossilIntList[f].fossilInt;
            if(fos->fossilIntList[f].prec != NOSUCH)
                fos->fossilIntList[index[f]].prec = index[fos->fossilIntList[f].prec];
            else
                fos->fossilIntList[index[f]].prec = NOSUCH;
        }
    fos->sizeFossil = size;
    free((void*)index);
}

void removeFossilIntNode(TypeFossilIntFeature *fos, int n) {
    if(fos == NULL)
        return;
    int *list, size=0, f;
    for(f=fos->fossilInt[n]; f != NOSUCH; f=fos->fossilIntList[f].prec)
        size++;
    if(size>0) {
        list = (int*) malloc((size+1)*sizeof(int));
        size = 0;
        for(f=fos->fossilInt[n]; f != NOSUCH; f=fos->fossilIntList[f].prec)
            list[size++] = f;
        list[size] = NOSUCH;
        removeFossilIntFossil(fos, list);
        free((void*)list);
    }
}

/*index[i] is the new index of former index i (NOSUCH if removed)*/
/*assume that no risk of ecrase a data*/
void reindexFossilInt(TypeFossilIntFeature *fos, int *index) {
    int n, *indexF, *indexT, ind;
    if(fos == NULL || index == NULL)
        return;
    indexF = (int*) malloc(fos->sizeFossil*sizeof(int));
    indexT = (int*) malloc(fos->sizeTime*sizeof(int));
    for(n=0; n<fos->sizeFossil; n++)
        indexF[n] = 1;
    for(n=0; n<fos->sizeTime; n++)
        indexT[n] = 1;
    ind = 0;
    for(n=0; n<fos->sizeNode; n++) {
        if(index[n] == NOSUCH) {
            if(fos->fossilInt[n] != NOSUCH) {
                int f;
                for(f=fos->fossilInt[n]; f!=NOSUCH; f = fos->fossilIntList[f].prec)
                    indexF[f] = NOSUCH;
                fos->fossilInt[n] = NOSUCH;
            }
            if(fos->endTime[n] != NOSUCH) {
                indexT[fos->endTime[n]] = NOSUCH;
                fos->endTime[n] = NOSUCH;
            }
        } else {
            if(fos->status[n] == unknownNodeStatus)
                warning("Warnin FoffilInt : unknown %d\n", n);
            fos->status[index[n]] = fos->status[n];
            fos->fossilInt[index[n]] = fos->fossilInt[n];
            fos->endTime[index[n]] = fos->endTime[n];
            ind++;
        }
    }
    fos->sizeNode = ind;
    ind=0;
    for(n=0; n<fos->sizeFossil; n++)
        if(indexF[n] != NOSUCH) {
            indexF[n] = ind++;
            fos->fossilIntList[indexF[n]] = fos->fossilIntList[n];
            if(fos->fossilIntList[indexF[n]].prec != NOSUCH)
                fos->fossilIntList[indexF[n]].prec = indexF[fos->fossilIntList[indexF[n]].prec];
        }
    fos->sizeFossil = ind;
    ind=0;
    for(n=0; n<fos->sizeTime; n++)
        if(indexT[n] != NOSUCH) {
            indexT[n] = ind++;
            fos->endTimeTable[indexT[n]] = fos->endTimeTable[n];
        }
    fos->sizeTime = ind;
    for(n=0; n<fos->sizeNode; n++) {
        if(fos->fossilInt[n] != NOSUCH)
            fos->fossilInt[n] = indexF[fos->fossilInt[n]];
        if(fos->endTime[n] != NOSUCH)
            fos->endTime[n] = indexT[fos->endTime[n]];
    }
    for(n=fos->sizeNode; n<fos->sizeBufNode; n++)
        fos->fossilInt[n] = NOSUCH;
    free((void*)indexF);
    free((void*)indexT);
}

void reallocFossilIntNode(TypeFossilIntFeature *fos, int size) {
    int i;
    if(fos == NULL)
        return;
    fos->fossilInt = (int*) realloc(fos->fossilInt, size*sizeof(int));
    for(i=fos->sizeBufNode; i<size; i++)
        fos->fossilInt[i] = NOSUCH;
    fos->status = (TypeNodeStatus*) realloc(fos->status, size*sizeof(TypeNodeStatus));
    for(i=fos->sizeBufNode; i<size; i++)
        fos->status[i] = noneNodeStatus;
    fos->endTime = (int*) realloc(fos->endTime, size*sizeof(int));
    for(i=fos->sizeBufNode; i<size; i++)
        fos->endTime[i] = NOSUCH;
    fos->sizeBufNode = size;
}

TypeIndexFossilIntTab *name2indexFossilIntTab(TypeNameFossilIntTab *tab, char **name, int size) {
    TypeLexiTree *dict;
    TypeIndexFossilIntTab *res;
    int i;
    res = (TypeIndexFossilIntTab *) malloc(sizeof(TypeIndexFossilIntTab));
    res->size = 0;
    res->index = (int*) malloc(tab->size*sizeof(int));
    res->fossilIntTab = (TypeFossilIntTab*) malloc(tab->size*sizeof(TypeFossilIntTab));
    dict = getDictNameTab(name, size);
    for(i=0; i<tab->size; i++) {
        int n;
        if((n = findWordLexi(tab->name[i], dict))>=0) {
            res->index[res->size] = n;
            res->fossilIntTab[res->size] = tab->fossilIntTab[i];
            res->size++;
        }
    }
    res->index = (int*) realloc((void*)res->index, res->size*sizeof(int));
    res->fossilIntTab = (TypeFossilIntTab*) realloc((void*)res->fossilIntTab, res->size*sizeof(double));
    freeLexiTree(dict);
    return res;
}



TypeFossilIntFeature *getFossilIntFeature(FILE *fi, char **name, int size) {
	TypeLexiTree *dict;
	TypeNameFossilIntTab *list;
	int i, n, f, *tmp;
	size_t *index;
	TypeFossilIntFeature *fos;
	dict = getDictNameTab(name, size);
	list = readFossilIntTab(fi);
	fos = (TypeFossilIntFeature *) malloc(sizeof(TypeFossilIntFeature));
	fos->endTime = NULL;
	fos->sizeFossil = 0;
	fos->sizeBufFossil = INC_FOSSIL_ITEM;
	fos->sizeBufNode = size;
	fos->sizeNode = size;
	fos->fossilIntList = (TypeFossilIntList*) malloc(fos->sizeBufFossil*sizeof(TypeFossilIntList));
	fos->fossilInt = (int*) malloc(size*sizeof(int));
	fos->status = (TypeNodeStatus*) malloc(size*sizeof(TypeNodeStatus));
	fos->endTimeTable = (TypeTimeInterval*)malloc(size*sizeof(TypeTimeInterval));
	for(i=0; i<size; i++) {
		fos->fossilInt[i] = NOSUCH;
		fos->status[i] = noneNodeStatus;
	}
	for(i=0; i<list->size; i++) {
		if((n = findWordLexi(list->name[i], dict))>=0) {
			TypeFossilIntTab ftab;
			ftab = list->fossilIntTab[i];
			if(list->stopTime[i] != NO_TIME) {
				int ind = 0;
				for(f=0; f<ftab.size; f++)
				if(ftab.fossilInt[f].inf>=list->stopTime[i] && ftab.fossilInt[f].sup>=list->stopTime[i])
					ftab.fossilInt[ind++] = ftab.fossilInt[f];
				ftab.size = ind;
				fos->status[n] = unknownNodeStatus;
				fos->endTimeTable[n].inf = list->stopTime[i];
				fos->endTimeTable[n].sup = list->stopTime[i];
			}
			if(ftab.size>0) {
				int prec;
				if(fos->sizeFossil+ftab.size >= fos->sizeBufFossil) {
					fos->sizeBufFossil += ftab.size+INC_FOSSIL_ITEM;
					fos->fossilIntList = (TypeFossilIntList*) realloc((void*)fos->fossilIntList, fos->sizeBufFossil*sizeof(TypeFossilIntList));
				}
				prec = -1;
				for(f=0; f<ftab.size; f++) {
					fos->fossilIntList[fos->sizeFossil].fossilInt = ftab.fossilInt[f];
					fos->fossilIntList[fos->sizeFossil].prec = prec;
					prec = fos->sizeFossil;
					fos->sizeFossil++;
				}
				fos->fossilInt[n] = fos->sizeFossil-1;
			}
		}
		free((void*) list->fossilIntTab[i].fossilInt);
	}
	for(i=0; i<list->size; i++)
		if(list->name[i] != NULL)
			free(list->name[i]);
	free(list->name);
	free(list->stopTime);
	free(list->fossilIntTab);
	free(list);
	freeLexiTree(dict);
	index = qsortindex(fos->fossilIntList, fos->sizeFossil, sizeof(TypeFossilIntList), compareFossilInt);
	for(n=0; n<size; n++)
		if(fos->fossilInt[n]>=0)
			fos->fossilInt[n] = index[fos->fossilInt[n]];
	for(f=0; f<fos->sizeFossil; f++)
		if(fos->fossilIntList[f].prec>=0)
			fos->fossilIntList[f].prec = index[fos->fossilIntList[f].prec];
	free((void*) index);
	tmp = (int*) malloc((fos->sizeFossil+1)*sizeof(int));
	for(n=0; n<size; n++) {
		int ind = 0;
		for(f=fos->fossilInt[n]; f>=0; f=fos->fossilIntList[f].prec)
			tmp[ind++] = f;
		if(ind>0) {
			qsort(tmp, ind, sizeof(int), compareInt);
			fos->fossilIntList[tmp[0]].prec = -1;
			for(f=1; f<ind; f++)
				fos->fossilIntList[tmp[f]].prec = tmp[f-1];
			fos->fossilInt[n] = tmp[ind-1];
		}
	}
	free((void*)tmp);
	return fos;
}

/*print fossilInt*/
void fprintFossilIntFeature(FILE *f, TypeFossilIntFeature *feat, char **name, int size) {
    int i;
    for(i=0; i<size; i++) {
        if(feat->fossilInt[i]>=0) {
            if(name != NULL && name[i] != NULL)
                fprintf(f, "%s", name[i]);
            else
                fprintf(f, "---");
            fprintFossilIntList(f, feat->fossilInt[i], feat->fossilIntList);
            fprintf(f, "\n");
        }
    }
}

/*print fossilInt*/
void fprintFossilInt(FILE *f, TypeTimeInterval fos) {
    if(fos.inf == fos.sup)
        fprintf(f, "%lf", fos.inf);
    else
        fprintf(f, "[%lf;%lf]", fos.inf, fos.sup);
}

/*print fossilInt table*/
void fprintFossilIntList(FILE *fo, int e, TypeFossilIntList *list) {
    int f;
    for(f=e; f>=0; f=list[f].prec) {
        fprintf(fo, "\t");
        fprintFossilInt(fo, list[f].fossilInt);
    }
}

/*print fossilInt table*/
void fprintFossilIntTabNHX(FILE *f, TypeFossilIntTab tab) {
    int i;
    if(!tab.size)
        return;
    for(i=0; i<tab.size; i++) {
        fprintf(f, ":FOS=");
        fprintFossilInt(f, tab.fossilInt[i]);
    }
}

/*print fossilInt table*/
void fprintFossilIntListNHX(FILE *fo, int e, TypeFossilIntList *list) {
    int f;
    for(f=e; f!=NOSUCH; f = list[f].prec) {
        fprintf(fo, ":FOS=");
        fprintFossilInt(fo, list[f].fossilInt);
    }
}


/*print fossilInt table*/
int sprintFossilIntListNHX(char *s, int e, TypeFossilIntList *list) {
    int f, ind = 0;
    for(f=e; f!=NOSUCH; f = list[f].prec) {
        ind += sprintf(s+ind, " :%s=", TAG_FOSSIL);
        ind += sprintInterval(s+ind, list[f].fossilInt);
    }
    return ind;
}

#define STRING_INTER_LENGTH 60
void fillComment(TypeTree *tree, TypeFossilIntFeature *fos) {
    if(tree == NULL || fos == NULL || tree->size == 0)
        return;
    if(tree->comment == NULL)
        tree->comment = (char**) malloc(tree->size*sizeof(char*));
    int n;
    for(n=0; n<tree->root; n++) {
        int f, nb = 0, length;
        switch(fos->status[n]) {
            case noneNodeStatus:
                nb = 0;
            break;
            case contempNodeStatus:
            case extinctNodeStatus:
                nb = 1;
            break;
            case unknownNodeStatus:
                nb = 3;
            break;
            default:
                nb = 0;
        }
        for(f=fos->fossilInt[n]; f!=NOSUCH; f = fos->fossilIntList[f].prec)
            nb++;
        if(tree->comment[n] != NULL)
            free((void*)tree->comment[n]);
        if(nb>0) {
            tree->comment[n] = (char*) malloc(nb*STRING_INTER_LENGTH*sizeof(char)+1);
             length = sprintFossilIntListNHX(tree->comment[n], fos->fossilInt[n], fos->fossilIntList);
            if(fos->status[n]==contempNodeStatus||fos->status[n]==extinctNodeStatus||fos->status[n]==unknownNodeStatus)
                length += sprintf(tree->comment[n]+length," :%s=%d", TAG_STATUS, (int) fos->status[n]);
            if(fos->status[n] == unknownNodeStatus) {
                length += sprintf(tree->comment[n]+length," :%s=", TAG_TIME);
                length += sprintInterval(tree->comment[n]+length, fos->endTimeTable[fos->endTime[n]]);
            }
            tree->comment[n][length] = '\0';
            tree->comment[n] = (char*) realloc(tree->comment[n], (length+1)*sizeof(char));
        } else
            tree->comment[n] = NULL;
    }
    {
        int f, nb = 0, length = 0;
        switch(fos->status[tree->root]) {
            case noneNodeStatus:
                nb = 0;
            break;
            case contempNodeStatus:
            case extinctNodeStatus:
                nb = 1;
            break;
            case unknownNodeStatus:
                nb = 3;
            break;
            default:
                nb = 0;
        }
        for(f=fos->fossilInt[tree->root]; f!=NOSUCH; f = fos->fossilIntList[f].prec)
            nb++;
        if(tree->comment[tree->root] != NULL)
            free((void*)tree->comment[tree->root]);
        if(tree->minTimeInt.inf!=NO_TIME && tree->minTimeInt.sup!=NO_TIME)
            nb++;
        if(tree->maxTimeInt.inf!=NO_TIME && tree->maxTimeInt.sup!=NO_TIME)
            nb++;
        if(nb>0) {
            tree->comment[tree->root] = (char*) malloc((nb)*STRING_INTER_LENGTH*sizeof(char)+1);
            if(tree->minTimeInt.inf!=NO_TIME && tree->minTimeInt.sup!=NO_TIME) {
                length += sprintf(tree->comment[tree->root]+length, " :ORI=");
                length += sprintInterval(tree->comment[tree->root]+length, tree->minTimeInt);
            }
            if(tree->maxTimeInt.inf!=NO_TIME && tree->maxTimeInt.sup!=NO_TIME) {
                length += sprintf(tree->comment[tree->root]+length, " :END=");
                length += sprintInterval(tree->comment[tree->root]+length, tree->maxTimeInt);
            }
            length += sprintFossilIntListNHX(tree->comment[tree->root]+length, fos->fossilInt[tree->root], fos->fossilIntList);
            if(fos->status[n]==contempNodeStatus||fos->status[tree->root]==extinctNodeStatus||fos->status[tree->root]==unknownNodeStatus)
                length += sprintf(tree->comment[tree->root]+length," :%s=%d", TAG_STATUS, (int) fos->status[tree->root]);
            if(fos->status[tree->root] == unknownNodeStatus) {
                length += sprintf(tree->comment[tree->root]+length," :%s=", TAG_TIME);
                length += sprintInterval(tree->comment[tree->root]+length, fos->endTimeTable[fos->endTime[tree->root]]);
            }
            tree->comment[tree->root][length] = '\0';
            tree->comment[tree->root] = (char*) realloc(tree->comment[tree->root], (length+1)*sizeof(char));
        } else
            tree->comment[tree->root] = NULL;
    }
    for(n=tree->root+1; n<tree->size; n++) {
        int f, nb = 0, length;
        switch(fos->status[n]) {
            case noneNodeStatus:
                nb = 0;
            break;
            case contempNodeStatus:
            case extinctNodeStatus:
                nb = 1;
            break;
            case unknownNodeStatus:
                nb = 3;
            break;
            default:
                nb = 0;
        }
        for(f=fos->fossilInt[n]; f!=NOSUCH; f = fos->fossilIntList[f].prec)
            nb++;
        if(tree->comment[n] != NULL)
            free((void*)tree->comment[n]);
        if(nb>0) {
            tree->comment[n] = (char*) malloc(nb*STRING_INTER_LENGTH*sizeof(char)+1);
             length = sprintFossilIntListNHX(tree->comment[n], fos->fossilInt[n], fos->fossilIntList);
            if(fos->status[n]==contempNodeStatus||fos->status[n]==extinctNodeStatus||fos->status[n]==unknownNodeStatus)
            length += sprintf(tree->comment[n]+length," :%s=%d", TAG_STATUS, (int) fos->status[n]);
            if(fos->status[n] == unknownNodeStatus) {
                length += sprintf(tree->comment[n]+length," :%s=", TAG_TIME);
//                length += sprintInterval(tree->comment[n]+length, fos->endTimeTable[fos->endTime[n]]);
                 length += sprintf(tree->comment[n]+length,"%lf", tree->time[n]);
           }
            tree->comment[n][length] = '\0';
            tree->comment[n] = (char*) realloc(tree->comment[n], (length+1)*sizeof(char));
        } else
            tree->comment[n] = NULL;
    }
}

/*Save tree with fossil intervals as comments*/
void fprintTreeFossilInt(FILE *f, TypeTree *tree, TypeFossilIntFeature *fos) {
    fprintSubTreeFossilInt(f, tree->root, tree, fos);
}

/*Save subtree at node with fossil intervals as comments*/
void fprintSubTreeFossilInt(FILE *f, int node, TypeTree *tree, TypeFossilIntFeature *fos) {
    char **comment_saved, **name_saved;
    int n;
    name_saved = tree->name;
    if(tree->name == NULL)
        tree->name = nameLeaves("leaf", tree);
    comment_saved = tree->comment;
    tree->comment = (char**) malloc(tree->sizeBuf*sizeof(char*));
    for(n=0; n<tree->sizeBuf; n++)
        tree->comment[n] = NULL;
    fillComment(tree, fos);
    fprintSubtreeNewick(f, node, tree);
    for(n=0; n<tree->sizeBuf; n++)
        if(tree->comment[n] != NULL)
            free((void*)tree->comment[n]);
    free((void*)tree->comment);
    tree->comment = comment_saved;
    if(name_saved == NULL) {
        for(n=0; n<tree->sizeBuf; n++)
            if(tree->name[n] != NULL)
                free((void*)tree->name[n]);
        free((void*)tree->name);
    }
    tree->name = name_saved;
}



void fillBoundsFossilInt(int n, double tmin, double tmax, TypeTree *tree,  TypeFossilIntFeature *fos, double *min, double *max, int *dmax) {
	int c;
	if(tree->time[n] != NO_TIME) {
		min[n] = tree->time[n];
	} else {
		if(fos && fos->fossilInt[n] != NOSUCH) {
			int f;
			min[n] = fos->fossilIntList[fos->fossilInt[n]].fossilInt.sup;
			for(f=fos->fossilIntList[fos->fossilInt[n]].prec; f!=NOSUCH; f=fos->fossilIntList[f].prec)
				if(fos->fossilIntList[f].fossilInt.sup>min[n])
					min[n] = fos->fossilIntList[f].fossilInt.sup;
		} else
			min[n] = tmin;
	}
	for(c=tree->node[n].child; c>=0; c=tree->node[c].sibling)
		fillBoundsFossilInt(c, min[n], tmax, tree,  fos, min, max, dmax);
	if(tree->time[n] != NO_TIME) {
		max[n] = tree->time[n];
		dmax[n] = 0;
	} else {
		if(tree->node[n].child == NOSUCH) {
			max[n] = tmax;
			dmax[n] = 0;
		} else {	
			max[n] = tmax+1;
			dmax[n] = 0;
			for(c=tree->node[n].child; c>=0; c=tree->node[c].sibling) {
				if(fos && fos->fossilInt[c] != -1) {
					int f;
					for(f=fos->fossilInt[c]; f>=0; f=fos->fossilIntList[f].prec)
						if(fos->fossilIntList[f].fossilInt.inf<max[n]) {
							max[n] = fos->fossilIntList[f].fossilInt.inf;
							dmax[n] = 0;
						}
				} else {
					if((max[c])<(max[n])) {
						max[n] = max[c];
						dmax[n] = dmax[c]+1;
					}
				}
			}
		}
	}
}

void fillUnknownTimesFossilInt(double tmin, double tmax, TypeTree *tree,  TypeFossilIntFeature *fos) {
	int *dmax;
	double *min, *max;
	
	min = (double*) malloc(tree->size*sizeof(double));
	max = (double*) malloc(tree->size*sizeof(double));
	dmax = (int*) malloc(tree->size*sizeof(int));
	fillBoundsFossilInt(tree->root, tmin, tmax, tree,  (TypeFossilIntFeature*) fos, min, max, dmax);
	fillTime(tree->root, tmin, tree, min, max, dmax);
	free((void*)min);
	free((void*)max);
	free((void*)dmax);
}
