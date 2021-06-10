#ifndef TreeF
#define TreeF
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include "Utils.h"




#define INC_SIZE 100
#define INC_SIZE_TREE 100
//#define MAX_CURRENT INT_MAX
#define MAX_CURRENT 50000
//#define MAX(x,y) ((x)>(y)?(x):(y))
//#define MIN(x,y) ((x)<(y)?(x):(y))
#define UNKNOWN -1
#define NOSUCH -1
#define NO_TIME -DBL_MAX
#define CO_TIME DBL_MAX

typedef enum DISPLAY_NAME {
	display_none=0,
	display_name,
	display_index,
	display_both,
	display_time_none,
	display_time_name,
	display_time_index,
	display_time_both
} TypeDisplayName;

typedef struct TIME_INTERVAL {
    double inf, sup;
} TypeTimeInterval;

typedef struct NODE {
    int child, sibling;
} TypeNode;

typedef struct TIME_TAB {
    char **name1, **name2;
    double *time;
    int size;
} TypeTimeTab;

typedef void* (*TypeGetFeature)(int, void*);
typedef void (*TypeSetFeature)(int, void*, void*);
typedef void* (*TypeAllocFeature)(int);
typedef void (*TypeReallocFeature)(int, void*);
typedef void (*TypeFreeFeature)(int, void*);

typedef struct FEATURE_FUNCT {
    TypeReallocFeature realloc;
    TypeFreeFeature free;
    TypeSetFeature set;
} TypeFeatureFunc;

typedef struct TREE_FEATURE {
    void *data;
    TypeAllocFeature alloc;
    TypeSetFeature set;
    TypeGetFeature get;
    TypeFreeFeature free;
} TypeTreeFrature;

/*
typedef struct STANDARD_FEATURE {
    char **name, **comment;
    int size, sizeBuf;
    void *info;
} TypeStandardFeature;
*/

typedef struct TREE {
    TypeNode *node;
    int root, size, sizeBuf, *parent;
    double *time, maxTime, minTime;
    char **name, **comment;
    void *info;
    TypeTimeInterval minTimeInt, maxTimeInt;
} TypeTree;


#ifdef __cplusplus
extern "C" {
namespace Tree {
#endif
/*print tree in newick format*/
void sprintTreeNewick(char *f, TypeTree *tree);double getMinTimeFromNode(int n, TypeTree *tree);
int getClade(char **listClad, TypeTree *tree);
int getNodeFromName(char *name, TypeTree *tree);
void fprintLeavesFromNode(FILE *f, int n, TypeTree *tree);
double getMaximumLeafTime(TypeTree *tree);
double getMaximumTime(TypeTree *tree);
void fillTime(int n, double tanc, TypeTree *tree, double *min, double *max, int *dmax);
void fillBounds(int n, double tmin, double tmax, TypeTree *tree, double *min, double *max, int *dmax);
void fillUnknownTimes(double tmin, double tmax, TypeTree *tree);
/*set the table of parents*/
/*remove the subtree pending from n*/
int *removeSubtreeReturnIndex(int n, TypeTree *tree);
void setParent(TypeTree *tree);
TypeTree *readTreeFromName(char *filename);
/*remove the subtree pending from n*/
void removeSubtree(int n, TypeTree *tree);
/*add a new leaf from branch/node n with name n*/
int addLeaf(int n, char *name, TypeTree *tree);
/*extend the tree buffer*/
void reallocTreeFeature(TypeTree *tree);
/*make tree to be binary*/
void toBinary(TypeTree *tree);
/*make node m be child of n*/
void transfer(int m, int n, TypeTree *tree);
/*return the root of tree*/
int getRoot(TypeTree *tree);
/*return the numbre of childre of node n*/
int getNumberChildren(int n, TypeTree *tree);
/*returns the number of leaves of the n subtree of "tree"*/
int countSubLeaves(int n, TypeTree *tree);
/*name (numerote) leaves of tree*/
char **nameBoth(char *prefixIntern, char *prefixLeaf, TypeTree *tree);
/*name (numerote) leaves of tree*/
char **nameLeaves(char *prefix, TypeTree *tree);
/*return a standard feature*/
//TypeStandardFeature *getBasicStandardFeature(TypeTree *tree);
/*turn branch lengthes to absolute time*/
void bltoabsTime(TypeTree *tree);
/*turn branch lengthes to absolute time*/
void abstoblTime(TypeTree *tree);
/*get the smallest time of the tree - makes sense only with absolute times*/
double getMinTime(TypeTree *tree);
/*get the greatest time of the tree - makes sense only with absolute times*/
double getMaxTime(TypeTree *tree);
/*shift all the times/branch from beg*/
void offsetTime(double beg, TypeTree *tree);
/*return the least common ancestor of n and m*/
int getLCA(int n, int m, TypeTree *tree);
/*return 1 if m descends from n*/
int isDescendant(int m, int n, TypeTree *tree);
/*return the table of parents*/
int *getParent(TypeTree *tree);
/*get the status from comments*/
char *getSpecy(char *str);
/*return name tab containing species found in comments*/
char **getNamestoSpecies(char **comment, TypeTree *tree);
/*return name tab containing species found in comments for leaves only*/
char **getNamestoSpeciesLeaves(char **comment, TypeTree *tree);
/*returns the number of contemporary lineages (i.e. living at maxTime) of "tree"*/
int countContemp(TypeTree *tree);
/*returns the number of leaves of "tree"*/
int countLeaves(TypeTree *tree);
/*fully duplicate "tree"*/
TypeTree *cpyTree(TypeTree *tree);
/*create and allocate a standard feature of sizeBuf*/
void *newStandardFeature(int sizeBuf);
/*reallocate a standard feature of sizeBuf*/
void reallocStandardFeature(int sizeBuf, void *std);
/*allocate and initialize a new_feat tree*/
TypeTree *newTree(int sizeBuf);
/*desallocate standard features*/
void freeStandardFeature(int size, void *std);
void setStandardFeature(int n, void *dest, int m, void *src);
/*desallocate tree*/
void freeTree(TypeTree *tree);
/*initialize fields of node to standard values*/
void initNode(TypeNode *node);
/*print tree in debug mode*/
void fprintTreeX(FILE *f, TypeTree *tree);
/*replace ' ' by underscores*/
void fixNameUnderscore(char **name, int size);
/*print ident, time and comment of node n*/
void fprintIdentTimeComment(FILE *f, int n, TypeTree *tree);
/*print tree in newick format*/
void fprintSubtreeNewick(FILE *f, int n, TypeTree *tree);
/*print tree in newick format*/
void fprintTreeNewick(FILE *f, TypeTree *tree);
/*print node in newick format*/
void fprintNodeNewick(FILE *f, int n, TypeTree *tree);
/*read tree in newick format*/
TypeTree *readTree(FILE *f);
/*read a node in newick format*/
int readNode(FILE *f, TypeTree *tree);
/*read a leaf in newick format*/
int readLeaf(FILE *f, TypeTree *tree);
/*read the name*/
char *readName(FILE *f);
/*read the comment*/
char *readComment(FILE *f);
/*read the time or branch length*/
double readTime(FILE *f);
/*reorder all children in a canonical (provided at least leaves are named)*/
void reorderTree(char **name, TypeTree *tree);
/*reorder all children in a canonical (provided at least leaves are named)*/
TypeTree *reorderTreeCpy(char **name, TypeTree *tree);
/*return a tree where leaves appear in lexicographic order*/
TypeLexiTree *getDictTree(char **name, int size);
/*print tree in newick format*/
void fprintTreePst(FILE *f, TypeTree *tree);
/*print node in newick format*/
void fprintNodePst(FILE *f, int n, TypeTree *tree);
/*rename internal nodes*/
void nameInternalNodes(char **name, TypeTree *tree);
TypeTimeTab *readTimeTab(FILE *f);
/*read a table <leaf1>\t<leaf2>\t<time> and set the time of the lca of <leaf1> and <leaf2> to <time>*/
void setTime(FILE *fi, char **name, TypeTree *tree);
/*prune "tree" to that can be observed from contemporary lineages only*/
TypeTree *pruneContemp(TypeTree *tree);
TypeTree *fixBinary(TypeTree *tree);
extern TypeFeatureFunc stdFeatureFunc;
/*print tree in text format*/
void printTreeDebug(FILE *f, int n, TypeTree *tree, char **name);
/*reallocate tree*/
void reallocTree(int sizeBuf, TypeTree *tree);
/*allocate a string tab and set its entries to NULL*/
char **newStringTab(int size);
TypeTree **readTrees(FILE *f);
void reorderTreeSize(TypeTree *tree);
void reindexTree(TypeTree *tree, int *index);
void fillIndexTree(int n, int *cur, TypeTree *tree, int *index);
/*make the set of trees 'tree' such that each node name has the same index in each of the tree, the nodes with names have the first indexes */
void fixTreeSet(TypeTree **tree);
#ifdef __cplusplus
}
}
#endif

#endif
