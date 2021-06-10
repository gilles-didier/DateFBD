#ifndef FossilIntF
#define FossilIntF

#include <stdlib.h>
#include <stdio.h>
#include "Utils.h"
#include "Tree.h"
#include "Fossil.h"




typedef struct FOSSIL_INT_LIST {
    TypeTimeInterval fossilInt;
    int prec;
} TypeFossilIntList;

typedef struct FOSSIL_INT_TAB {
    TypeTimeInterval *fossilInt;
    int size;
} TypeFossilIntTab;

typedef struct INDEX_FOSSIL_INT_TAB {
    int size, *index;
    TypeFossilIntTab *fossilIntTab;
} TypeIndexFossilIntTab;

typedef struct NAME_FOSSIL_INT_TAB {
    char **name;
    TypeFossilIntTab *fossilIntTab;
    double *stopTime;
    int size;
} TypeNameFossilIntTab;

typedef struct FOSSIL_INT_FEATURE {
    TypeFossilIntList *fossilIntList;
    TypeTimeInterval *endTimeTable;
    int *fossilInt, sizeFossil, sizeBufFossil, sizeTime, sizeBufTime, sizeNode, sizeBufNode, *endTime;
    TypeNodeStatus *status;
} TypeFossilIntFeature;

typedef struct TMP_FOSSIL_INT_ITEM {
    char *name;
    TypeFossilIntTab fossilIntTab;
} TypeFossilIntItem;

#ifdef __cplusplus
extern "C" {
#endif
double getFossilIntMaxInfTime(TypeFossilIntFeature* feat);
double getFirstFossilIntTimeFromNode(int n, TypeTree *tree, TypeFossilIntFeature* feat);
double getMidFossilIntTimeFromNode(int n, TypeTree *tree, TypeFossilIntFeature* feat);
/*index[i] is the new index of former index i (NOSUCH if removed)*/
/*assume that no risk of ecrase a data*/
void reindexFossilInt(TypeFossilIntFeature *fos, int *index);
/*fix status of leaves*/
void fixStatus(TypeTree *tree, TypeFossilIntFeature *feat);
void fillComment(TypeTree *tree, TypeFossilIntFeature *fos);
/*print fossilInt table*/
int sprintFossilIntListNHX(char *s, int e, TypeFossilIntList *list);
/*from fos to FossilInt fos*/
TypeFossilIntFeature *fosToFossilInt(TypeTree *tree);
void removeFossilIntNode(TypeFossilIntFeature *fos, int n);
void removeFossilIntFossil(TypeFossilIntFeature *fos, int* list);
void reallocFossilIntNode(TypeFossilIntFeature *fos, int size);
/*allocate a new fossil int feature with n entries/nodes*/
TypeFossilIntFeature *newFossilIntFeature(int n);
/*add a new fossil to n*/
void addFossilInt(TypeFossilIntFeature *fos, int n, double inf, double sup);
TypeFossilIntFeature *cpyFossilIntFeature(TypeFossilIntFeature *feat, int n);
void fixTree(TypeTree *tree, TypeFossilIntFeature *fos);
TypeTimeInterval toFossilInt(char *s);
TypeFossilIntTab getFossilInt(char *str);
TypeFossilIntFeature *treeToFossilInt(TypeTree *tree);
TypeNameFossilIntTab *readFossilIntTab(FILE *f);
void freeNameFossilIntTab(TypeNameFossilIntTab *list);
void freeIndexFossilIntTab(TypeIndexFossilIntTab *list);
TypeIndexFossilIntTab *name2indexFossilIntTab(TypeNameFossilIntTab *tab, char **name, int size);
TypeFossilIntFeature *getFossilIntFeature(FILE *f, char **name, int size);
void fprintFossilInt(FILE *f, TypeTimeInterval fos);
void fprintFossilIntFeature(FILE *f, TypeFossilIntFeature *feat, char **name, int size);
void fprintFossilIntList(FILE *fo, int e, TypeFossilIntList *list);
TypeFossilFeature *sampleFossilInt(TypeFossilIntFeature* feat, int size);
double getFirstFossilIntTime(TypeFossilIntFeature* feat);
double getMidFossilIntTime(TypeFossilIntFeature* feat);
double getMinFossilIntTime(TypeFossilIntFeature* feat);
double getMaxFossilIntTime(TypeFossilIntFeature* feat);
double getMinFossilIntTimeFromNode(int n, TypeTree *tree, TypeFossilIntFeature* feat);
double getMaxFossilIntTimeToNode(int n, TypeTree *tree, TypeFossilIntFeature* feat);
void negateFossilInt(TypeFossilIntFeature* feat);
int compareFossilInt(const void* a, const void* b);
int compareFossilIntList(const void* a, const void* b);
void freeFossilIntFeature(TypeFossilIntFeature *fos);
void fprintTreeFossilInt(FILE *f, TypeTree *tree, TypeFossilIntFeature *fos);
/*Save subtree at node with fossil intervals as comments*/
void fprintSubTreeFossilInt(FILE *f, int node, TypeTree *tree, TypeFossilIntFeature *fos);
void fillUnknownTimesFossilInt(double tmin, double tmax, TypeTree *tree,  TypeFossilIntFeature *fos);
TypeFossilFeature *setTimeFossilInt(TypeFossilIntFeature* feat, double *time, int size);
#ifdef __cplusplus
}
#endif

#endif
