#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "Utils.h"
#define INC_BUFFER_UTILS 10

int DEBUGX=0;

static void fprintLexiTreeRec(FILE *f, char *tmp, int i, int n, TypeLexiTree *dict);

int isequal(char c) {
    return c == '=' || c ==':' || c == '~';
}

/*get the fossilInts from comments*/
char *nextTag(char *s, char **tag, char **value) {
    int i;

    for(i=0; s[i]!='\0' && s[i] != ':'; i++)
        ;
    if(s[i]!='\0') {
        i++;
        *tag = &(s[i]);
    } else {
        *tag = NULL;
        *value = NULL;
        return &(s[i]);
    }
    for(; s[i]!='\0' && !isspace(s[i]) && s[i] != '='; i++)
        ;
    if(s[i] != '=') {
        s[i++] = '\0';
        for(; s[i]!='\0' && isspace(s[i]) && s[i] != '='; i++)
            ;
        if(s[i] != '=') {
            *value = NULL;
            return &(s[i]);
        }
    } else
        s[i++] = '\0';
    for(; s[i]!='\0' && isspace(s[i]); i++)
        ;
    if(s[i]!='\0' && s[i]!=':')
        *value = &(s[i]);
    else {
        *value = NULL;
        return &(s[i]);
    }
   if(s[i] == '\'' || s[i] == '"') {
       for(; s[i]!='\0' && !(s[i] == '\'' || s[i] == '"'); i++)
           ;
       if(s[i] == '\'' || s[i] == '"')
           s[i++] = '\0';
       else
           *value = NULL;
       return &(s[i]);
    } else {
       for(; s[i]!='\0' && !isspace(s[i]); i++)
           ;
       if(s[i]!='\0')
            s[i++] = '\0';
       return &(s[i]);
   }
}



void numprint(char *buffer, double val) {
    sprintf(buffer, "%.2lf", val);
}

/*return a lexitree from tab name*/
TypeLexiTree *getDictNameTab(char **name, int size) {
    int i;
    TypeLexiTree *dict;
    dict = newLexiTree();
    for(i=0; i<size; i++)
        if(name && name[i]) {
            if(addWordLexi(name[i], i, dict)>=0)
                warning("Warning! duplicate identifier '%s'\n", name[i]);
        }
    return dict;
}

void freeListDouble(TypeListDouble *list) {
        free((void*)list->val);
        free((void*)list);
    }

void replaceChar(char *s, char a, char n) {
    int i;
    for(i=0; s[i] != '\0'; i++)
        if(s[i] == a)
            s[i] = n;
}

void removeSpaces(char *s) {
    int start, end, i, ind;
    for(start=0; s[start] != '\0' && issep(s[start]); start++);
    for(end=strlen(s)-1; end>0 && issep(s[end]); end--);
    ind = 0;
    for(i=start; i<=end; i++)
        s[ind++] = s[i];
    s[end+1] = '\0';
}

char *strdpl(char *src) {
    if(src != NULL) {
        char *res;
        res = (char*) malloc((strlen(src)+1)*sizeof(char));
        strcpy(res, src);
        return res;
    } else
        return NULL;
}

void fprintLexiTreeRec(FILE *f, char *tmp, int i, int n, TypeLexiTree *dict) {
    int c;
    for(c=dict->node[n].child; c>=0; c=dict->node[c].sibling) {
        if(dict->node[c].symbol == '\0') {
            tmp[i] = '\0';
            fprintf(f, "%s\t%d\n", tmp, dict->node[c].child);
        } else {
            tmp[i] = dict->node[c].symbol;
            fprintLexiTreeRec(f, tmp, i+1, c, dict);
        }
    }
}

void fprintLexiTree(FILE *f, TypeLexiTree *dict) {
    char tmp[MAX_DICT_LENGTH];
    fprintLexiTreeRec(f, tmp, 0, dict->root, dict);
}

int findWordLexi(char *w, TypeLexiTree *dict) {
    int i, n;
    n = dict->root;
    i = 0;
    while(i<MAX_DICT_LENGTH && w[i] != '\0') {
        int c;
        c = dict->node[n].child;
        while(c>=0 && dict->node[c].symbol<w[i])
            c = dict->node[c].sibling;
        if(c>=0 && dict->node[c].symbol == w[i]) {
            n = c;
            i++;
        } else
            return -1;
    }
    if(dict->node[n].child >= 0 && dict->node[dict->node[n].child].symbol == '\0')
        return dict->node[dict->node[n].child].child;
    else
        return -1;
}

/*add word w in dict. If w is already in, then it returns its index (stored in the child field of a leaf - labelled by '\0'), -1 otherwise*/
int addWordLexi(char *w, int index, TypeLexiTree *dict) {
    int i, n;
    n = dict->root;
    i = 0;
    while(i<MAX_DICT_LENGTH && w[i] != '\0') {
        int *prec, c;
        prec = &(dict->node[n].child);
        c = dict->node[n].child;
        while(c>=0 && dict->node[c].symbol<w[i]) {
            prec = &(dict->node[c].sibling);
            c = dict->node[c].sibling;
        }
        if(c>=0 && dict->node[c].symbol == w[i]) {
            n = c;
            i++;
        } else {
            *prec = dict->size;
            if(dict->size >= dict->sizeBuf) {
                dict->sizeBuf += INC_SIZE_DICT;
                dict->node = (TypeLexiNode*) realloc((void*)dict->node, dict->sizeBuf*sizeof(TypeLexiNode));
            }
            initLexiNode(w[i], &(dict->node[dict->size]));
            dict->node[dict->size].sibling = c;
            n = dict->size;
            dict->size++;
            i++;
            while(i<MAX_DICT_LENGTH && w[i] != '\0') {
                if(dict->size >= dict->sizeBuf) {
                    dict->sizeBuf += INC_SIZE_DICT;
                    dict->node = (TypeLexiNode*) realloc((void*)dict->node, dict->sizeBuf*sizeof(TypeLexiNode));
                }
                initLexiNode(w[i], &(dict->node[dict->size]));
                dict->node[n].child = dict->size;
                n = dict->node[n].child;
                dict->size++;
                i++;
            }
        }
    }
    if(dict->node[n].child >= 0 && dict->node[dict->node[n].child].symbol == '\0')
        return dict->node[dict->node[n].child].child;
    if(dict->size >= dict->sizeBuf) {
        dict->sizeBuf += INC_SIZE_DICT;
        dict->node = (TypeLexiNode*) realloc((void*)dict->node, dict->sizeBuf*sizeof(TypeLexiNode));
    }
    initLexiNode('\0', &(dict->node[dict->size]));
    dict->node[dict->size].sibling = dict->node[n].child;
    dict->node[dict->size].child = index;
    dict->node[n].child = dict->size;
    dict->size++;
    return -1;
}

void initLexiNode(char symbol, TypeLexiNode *n) {
    n->symbol = symbol;
    n->child = -1;
    n->sibling = -1;
}

TypeLexiTree *newLexiTree() {
    TypeLexiTree* dict;
    dict = (TypeLexiTree*) malloc(sizeof(TypeLexiTree));
    dict->sizeBuf = INC_SIZE_DICT;
    dict->node = (TypeLexiNode*) malloc(dict->sizeBuf*sizeof(TypeLexiNode));
    dict->root = 0;
    initLexiNode(0,&(dict->node[dict->root]));
    dict->size = 1;
    return dict;
}

void freeLexiTree(TypeLexiTree *dict) {
    if(dict->node != NULL)
        free((void*)dict->node);
    free((void*) dict);
}

void printIndex(FILE *f, TypeIndex *index) {
    int i;
    for(i=0; i<index->size; i++)
        fprintf(f, "%d\t%s\n", i, index->name[i]);
}


void initIndex(TypeIndex *index) {
    index->size = 0;
    index->buffer = INC_BUFFER_UTILS;
    index->name = (char **) malloc(index->buffer*sizeof(char*));
    index->dict = newDictNode('x');
}

int addIndex(char *name, TypeIndex *species) {
    int sizeTmp, index;
    if(strlen(name) == 0)
        return -1;
    sizeTmp = species->size;
    index = getIndexString(name, species->dict, &(species->size));
    if(species->size > sizeTmp) {
        if(sizeTmp >= species->buffer) {
            species->buffer += INC_BUFFER_UTILS;
            species->name = (char**) realloc(species->name, species->buffer*sizeof(char*));
        }
        species->name[index] = (char*) malloc((strlen(name)+1)*sizeof(char));
        strcpy(species->name[index], name);
    }
    return index;
}

TypeDictNode *newDictNode(char c){
    TypeDictNode *n;
    n = (TypeDictNode *) malloc(sizeof(TypeDictNode));
    n->symbol = c;
    n->child = NULL;
    n->sibling = NULL;
    n->index = -1;
    return n;
}

int getIndex(char *s, TypeDictNode *cur) {
    int i;

    for(i=0; s[i]!='\0'; i++) {
        TypeDictNode *tmp = cur->child;
        while(tmp != NULL && tmp->symbol < s[i])
            tmp = tmp->sibling;
        if(tmp == NULL || tmp->symbol > s[i])
            return -1;
        else
            cur = tmp;
    }
    return cur->index;
}

int getIndexString(char *s, TypeDictNode *cur, int *size) {
    int i;

    for(i=0; s[i]!='\0'; i++) {
        TypeDictNode *tmp;
        tmp = cur->child;
        if(tmp == NULL || tmp->symbol >= s[i]) {
            if(tmp == NULL || tmp->symbol > s[i]) {
                cur->child = newDictNode(s[i]);
                (cur->child)->sibling = tmp;
            }
            cur = cur->child;
        } else {
            while(tmp->sibling != NULL && (tmp->sibling)->symbol < s[i])
                tmp = tmp->sibling;
            if(tmp->sibling != NULL && (tmp->sibling)->symbol == s[i]) {
                cur = tmp->sibling;
            } else {
                TypeDictNode *retmp;
                retmp = tmp->sibling;
                tmp->sibling = newDictNode(s[i]);
                (tmp->sibling)->sibling = retmp;
                cur = tmp->sibling;
            }
        }
    }
    if(cur->index == -1)
        cur->index = (*size)++;
    return cur->index;
}


int IsSeparator(char c) {
    return (c == ' ' || c == '\t' || c == '\n' || c == '\r');
}
int IsItemSeparator(char c) {
    return (c == ' ' || c == '\t' || c == ';');
}
int IsLineSeparator(char c) {
    return (c == '\n' || c == '\r');
}

int issep(char c) {
    return (c == ' ' || c == '\t' || c == ';');
}


int isline(char c) {
    return (c == '\n' || c == '\r');
}

int issepline(char c) {
    return (c == ' ' || c == '\t' || c == ';' || c == '\n' || c == '\r');
}


void skipSeparator(FILE *f) {
    char c;
    while( (c = fgetc(f)) != EOF )
        if(!IsSeparator(c)) {
            ungetc(c, f);
            break;
        }
}

int readLine(FILE *f, char *buffer) {
    char c;
    int tot = 0;
    while((c=getc(f))!=EOF && IsLineSeparator(c));
    while(c!=EOF && !IsLineSeparator(c)) {
        buffer[tot++] = c;
        c=getc(f);
    }
    if(c!=EOF)
        buffer[tot++] = '\n';
    buffer[tot] = '\0';
    return tot;
}

int readItem(FILE *f, char *buffer) {
    char c;
    int tot = 0;
    while((c=getc(f))!=EOF && IsSeparator(c));
    while(c!=EOF && !IsSeparator(c)) {
        buffer[tot++] = c;
        c=getc(f);
    }
    buffer[tot++] = '\n';
    buffer[tot] = '\0';
    return tot;
}

char passLines(FILE *f, char c) {
    while(c != EOF && isspace(c))
        c = getc(f);
    while(c != EOF && (c == '%' || c == '#')) {
        while(c != EOF && !IsLineSeparator(c))
            c = getc(f);
        while(c != EOF && isspace(c))
            c = getc(f);
    }
    return c;
}


char nextStartLine(FILE *f, char c) {
    do {
        while(c != EOF && !IsLineSeparator(c))
            c = getc(f);
        while(c != EOF && isspace(c))
            c = getc(f);
    } while(c != EOF && (c == '%' || c == '#'));
    return c;
}

char nextStartItem(FILE *f, char c) {
        while(c != EOF && !IsItemSeparator(c)) {
            c = getc(f);
        }
        while(c != EOF && IsItemSeparator(c))
            c = getc(f);
    return c;
}

char passSpaces(FILE *f, char c) {
        while(c != EOF && IsItemSeparator(c))
            c = getc(f);
    return c;
}

char *truncFileName(char* name) {
    int i;
    for(i=strlen(name)-1; i>=0 && name[i] != '.'; i--);
    if(i>=0)
     name[i] = '\0';
    return name;
}

char *getExtension(char* name) {
    int i, length;

    length = strlen(name);
    for(i=length-1; i>0 && name[i] != '.'; i--);
    if(i>0)
        return &(name[i+1]);
    else
        return &(name[length]);
}

int tokenize(char *src, char **dest) {
    int indice = 0, pos = 0;
    do {
        while(src[pos] != '\0' && IsSeparator(src[pos]))
             pos++;
        if(src[pos] != '\0') {
            dest[indice++] = &(src[pos]);
        }
        while(src[pos] != '\0' && !IsSeparator(src[pos]))
             pos++;
        if(src[pos] != '\0') {
            src[pos++] = '\0';
        }
    } while(src[pos] != '\0');
    dest[indice] = NULL;
    return indice;
}

int find(char *src, char **dest, int size) {
    int i;
    for(i=0; i<size && strcmp(src, dest[i]) != 0; i++);
    return i;
}

void fixSpace(char *src) {
    int i;
    for(i=0; src[i] != '\0'; i++)
        if(IsSeparator(src[i]))
            src[i] = '_';
}

char skipLineSpaceComment(FILE *f, char c) {
    while(c != EOF && (isspace(c) || c == '#')) {
        while(c != EOF && isspace(c)) {
            c = getc(f);
        }
        if(c == '#') {
            do {
                c = getc(f);
            } while(c != EOF && !isline(c));
        }
    }
    return c;
}

char skipSep(FILE *f, char c) {
    while(c != EOF && issep(c)) {
        c = getc(f);
    }
    return c;
}

static int (*comparisonFunction)(const void *, const void*);

int compareInt(const void* a, const void* b) {
    if(*((int*)a)>*((int*)b))
        return 1;
    if(*((int*)a)<*((int*)b))
        return -1;
    return 0;
}

int compareIndirect(const void* a, const void* b) {
    const void *term1 = *((void**)a), *term2 = *((void**)b);
    return comparisonFunction(term1,term2);
}

/*sort the table base and return table index where entry i contains the index in the sorted table of the entry which was at i*/
size_t *qsortindex(void *base, size_t nitems, size_t size, int (*compar)(const void *, const void*)) {
    size_t *index, i;
    void **p, *data;
    if(size == 0 || nitems == 0)
        return NULL;
    index = (size_t*) malloc(nitems*sizeof(size_t*));
    p = (void**) malloc(nitems*sizeof(void*));
    data = (void*) malloc(nitems*size);
    memcpy(data, base, nitems*size);
    for(i=0; i<nitems; i++)
        p[i] = (void*)(((size_t)data)+i*size);
    comparisonFunction = compar;
    qsort(p, nitems, sizeof(void*), compareIndirect);
    for(i=0; i<nitems; i++) {
        index[(((size_t)p[i])-((size_t)data))/size] = i;
        memcpy((void*)(((size_t)base)+i*size), p[i], size);
    }
    free((void*)data);
    free((void*)p);
    return index;
}

/*return the ordered table of index of base*/
size_t *qsortTable(void *base, size_t nitems, size_t size, int (*compar)(const void *, const void*)) {
    size_t *index, i;
    void **p, *data;
    if(size == 0 || nitems == 0)
        return NULL;
    index = (size_t*) malloc(nitems*sizeof(size_t*));
    p = (void**) malloc(nitems*sizeof(void*));
    data = (void*) malloc(nitems*size);
    memcpy(data, base, nitems*size);
    for(i=0; i<nitems; i++)
        p[i] = (void*)(((size_t)data)+i*size);
    comparisonFunction = compar;
    qsort(p, nitems, sizeof(void*), compareIndirect);
    for(i=0; i<nitems; i++)
        index[i] = (((size_t)p[i])-((size_t)data))/size;
    free((void*)data);
    free((void*)p);
    return index;
}

int compareDouble(const void* a, const void* b) {
    if(*((double*)a)>*((double*)b))
        return 1;
    if(*((double*)a)<*((double*)b))
        return -1;
    return 0;
}

#define MAX_SIZE 100

void toInterval(char *s, double *inf, double *sup) {
    int i, ind;
    char tmp[MAX_SIZE+1];
    for(i=0; s[i] != '\0'; i++)
        if(s[i] == ',')
            s[i] = '.';
    for(i=0; s[i] != '\0' && isspace(s[i]); i++);
    if(s[i] == '(' || s[i] == '[') {
        i++;
        for(; s[i] != '\0' && isspace(s[i]); i++);
        ind = 0;
        for(; s[i] != '\0' && s[i] != ':' && s[i] != ';' && !isspace(s[i]) && ind<MAX_SIZE; i++)
            tmp[ind++] = s[i];
        tmp[ind] = '\0';
        if(ind>=MAX_SIZE)
            error("Error when reading an interval - number too long:\n%s...\n", tmp);
        *inf = atof(tmp);
        for(; s[i] != '\0' && isspace(s[i]); i++);
        if(s[i] != ':' && s[i] != ';')
            error("Missing ':' or ';' while reading an interval.");
        i++;
        for(; s[i] != '\0' && isspace(s[i]); i++);
        ind = 0;
        for(; s[i] != '\0' && s[i] != ')' && s[i] != ']' && !isspace(s[i]); i++)
            tmp[ind++] = s[i];
        tmp[ind] = '\0';
        if(ind>=MAX_SIZE)
            error("Error when reading an interval - number too long:\n%s...\n", tmp);
        *sup =	atof(tmp);
        for(; s[i] != '\0' && isspace(s[i]); i++);
        if(s[i] != ')' && s[i] != ']')
            error("Missing ')' or ']' while reading an interval.\n%s\n", s);
        if(*sup<*inf) {
            double tmp = *sup;
            *sup = *inf;
            *inf = tmp;
        }
    } else {
        *inf = atof(s);
        *sup = *inf;
    }
}
