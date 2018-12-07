#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <gsl/gsl_sf_gamma.h>
#include "Utils.h"
#include "TreeExtras.h"


//#define binomial(k,n) (gsl_sf_fact((unsigned int) n)/(gsl_sf_fact(((unsigned int) k))*gsl_sf_fact(((unsigned int) n-k))))



typedef struct TREE_RANKING {
    double rank;
    int leaf;
} TypeTreeRanking;

TypeTreeRanking logNumberRankingsRec(int n, TypeTree *tree);
TypeTreeRanking logNumberOrdersRec(int n, TypeTree *tree);
static double binomial(unsigned int k, unsigned int n);
static double logBinomial(unsigned int k, unsigned int n);

double binomial(unsigned int k, unsigned int n) {
    return exp(logBinomial(k,n));
}

double logBinomial(unsigned int k, unsigned int n) {
    return gsl_sf_lnfact(n)-(gsl_sf_lnfact(k)+gsl_sf_lnfact(n-k));
}


double numberRankings(TypeTree *tree) {
    TypeTreeRanking r = logNumberRankingsRec(tree->root, tree);
    return exp(r.rank);
}

double probTreeShape(TypeTree *tree) {
    return probSubTreeShape(tree->root, tree);
}

double probSubTreeShape(int n, TypeTree *tree) {
    return exp(logProbSubTreeShape(n, tree));
}

double logProbTreeShape(TypeTree *tree) {
    return logProbSubTreeShape(tree->root, tree);
}

double logProbSubTreeShape(int n, TypeTree *tree) {
    TypeTreeRanking  r = logNumberRankingsRec(n, tree);
    return r.rank-2*gsl_sf_lnfact(r.leaf-1)-log(r.leaf);
}

TypeTreeRanking logNumberRankingsRec(int n, TypeTree *tree) {
    TypeTreeRanking  rl, rr, res;
    if(tree->node[n].child == NOSUCH) {
        res.rank = 0.;
        res.leaf = 1;
        return res;
    }
    rl = logNumberRankingsRec(tree->node[n].child, tree);
    if(tree->node[tree->node[n].child].sibling == NOSUCH) {
        fprintf(stderr, "Warning: node with a single child\n");
        return rl;
    }
    rr = logNumberRankingsRec(tree->node[tree->node[n].child].sibling, tree);
    if(tree->node[tree->node[tree->node[n].child].sibling].sibling != NOSUCH)
        fprintf(stderr, "Warning: we don't deal with polytomies\n");
    res.leaf = rl.leaf+rr.leaf;
    res.rank = rl.rank+rr.rank+log(2)+logBinomial((unsigned int) (rl.leaf-1), (unsigned int) (res.leaf-2));
    return res;
}


double logNumberRankings(TypeTree *tree) {
    TypeTreeRanking r = logNumberRankingsRec(tree->root, tree);
    return r.rank;
}

TypeTreeRanking logNumberOrdersRec(int n, TypeTree *tree) {
    TypeTreeRanking  rl, rr, res;
    if(tree->node[n].child == NOSUCH) {
        res.rank = 0.;
        res.leaf = 1;
        return res;
    }
    rl = logNumberOrdersRec(tree->node[n].child, tree);
    if(tree->node[tree->node[n].child].sibling == NOSUCH) {
        fprintf(stderr, "Warning: node with a single child\n");
        return rl;
    }
    rr = logNumberOrdersRec(tree->node[tree->node[n].child].sibling, tree);
    if(tree->node[tree->node[tree->node[n].child].sibling].sibling != NOSUCH)
        fprintf(stderr, "Warning: we don't deal with polytomies\n");
    res.leaf = rl.leaf+rr.leaf;
    res.rank = rl.rank+rr.rank+logBinomial((unsigned int) (rl.leaf-1), (unsigned int) (res.leaf-2));
    return res;
}

double logNumberOrders(TypeTree *tree) {
    TypeTreeRanking r = logNumberOrdersRec(tree->root, tree);
    return r.rank;
}

double logSubNumberOrders(int n, TypeTree *tree) {
    TypeTreeRanking r = logNumberOrdersRec(n, tree);
    return r.rank;
}

double numberTrees(int n) {
	return exp(logNumberTrees(n));
}

double logNumberTrees(int n) {
    if(n < 2)
        return 0.;
    else
        return gsl_sf_lnfact(2*n-3)-gsl_sf_lnfact(n-2)-log(2.0)*((double) n-2.);
}



int compareEmpiricalItem(const void* a, const void* b) {
    if(((TypeEmpiricalDistributionItem*)a)->val>((TypeEmpiricalDistributionItem*)b)->val)
        return 1;
    if(((TypeEmpiricalDistributionItem*)a)->val<((TypeEmpiricalDistributionItem*)b)->val)
        return -1;
    return 0;
}


#define INC_BUFFER_DISTRIB 500
TypeEmpiricalDistributionItem **rankingNumberDistribution(int n) {
    int i;
    TypeEmpiricalDistributionItem **emp;
    emp = (TypeEmpiricalDistributionItem**) malloc(n*sizeof(TypeEmpiricalDistributionItem*));
    emp--;
    emp[1] = (TypeEmpiricalDistributionItem*) malloc(2*sizeof(TypeEmpiricalDistributionItem));
    emp[1][0].val = 0.;
    emp[1][0].eff = 1.;
    emp[1][1].val = sqrt(-1);
    for(i=2; i<=n; i++) {
        int j, ind = 0, size = INC_BUFFER_DISTRIB;
        emp[i] = (TypeEmpiricalDistributionItem*) malloc(size*sizeof(TypeEmpiricalDistributionItem));
        for(j=1; j<i; j++) {
            int k, l;
            for(k=0; !isnan(emp[j][k].val); k++)
                for(l=0; !isnan(emp[i-j][l].val); l++) {
                    if(ind >= size) {
                        size += INC_BUFFER_DISTRIB;
                        emp[i] = (TypeEmpiricalDistributionItem*) realloc((void*) emp[i], size*sizeof(TypeEmpiricalDistributionItem));
                    }
                    emp[i][ind].val = emp[j][k].val+emp[i-j][l].val+log(2)+logBinomial((unsigned int) (j-1), (unsigned int) (i-2));
                    emp[i][ind].eff = emp[j][k].eff*emp[i-j][l].eff*binomial((unsigned int) j, (unsigned int) i)/2.;
                    ind++;
                 }
         }

double suma=0., sumb=0.;
for(j=0; j<ind; j++) {
	suma+=emp[i][j].eff;
	sumb+=emp[i][j].eff*exp(emp[i][j].val);
}
printf("i %d tot eff %.2lE/%.2lE %.2lE/%.2lE\n", i, suma, numberTrees(i), sumb, exp(2*gsl_sf_lnfact(i-1)+log(i)));
        size = ind;
        qsort(emp[i], size, sizeof(TypeEmpiricalDistributionItem),compareEmpiricalItem);
       ind=0;
        for(j=1; j<size; j++) {
            if(emp[i][j].val == emp[i][ind].val)
                emp[i][ind].eff += emp[i][j].eff;
            else
                emp[i][++ind] = emp[i][j];
        }
        emp[i][++ind].val = sqrt(-1);
       emp[i] = (TypeEmpiricalDistributionItem*) realloc((void*) emp[i], (ind+1)*sizeof(TypeEmpiricalDistributionItem));
    }
    return emp+1;
}

double *getThresholds(double prop, TypeEmpiricalDistributionItem **emp, int n) {
	double *threshold;
	int i;
	threshold = (double*) malloc(n*sizeof(double));
	for(i=0; i<n; i++) {
		double sumE = 0., sumV = 0., cumE = 0., cumV = 0.;
		int j, size = 0;
		for(j=0; !isnan(emp[i][j].val); j++) {
			sumE += emp[i][j].eff;
			sumV += emp[i][j].eff*exp(emp[i][j].val);
			size++;
		}
		for(j=size-1; j>=0 && cumV/sumV < prop; j--) {
			cumE += emp[i][j].eff;
			cumV += emp[i][j].eff*exp(emp[i][j].val);
		}
		threshold[i] = cumE/sumE;
	}
	return threshold;
}

double *getThresholdsBis(double prop, TypeEmpiricalDistributionItem **emp, int n) {
	double *threshold;
	int i;
	threshold = (double*) malloc(n*sizeof(double));
	for(i=0; i<n; i++) {
		double sumE = 0., sumV = 0., cumE = 0., cumV = 0.;
		int j, size = 0;
		for(j=0; !isnan(emp[i][j].val); j++) {
			sumE += emp[i][j].eff;
			sumV += emp[i][j].eff*exp(emp[i][j].val);
			size++;
		}
		for(j=size-1; j>=0 && cumE/sumE < prop; j--) {
			cumE += emp[i][j].eff;
			cumV += emp[i][j].eff*exp(emp[i][j].val);
		}
		threshold[i] = cumV/sumV;
	}
	return threshold;
}

void toDistribution(TypeEmpiricalDistributionItem **emp, int n) {
    int i;
    for(i=0; i<n; i++) {
        int j;
/*        double sum = 0;
        for(j=0; !isnan(emp[i][j].val); j++)
            sum += emp[i][j].eff;
*/
      double sum = 0;
        for(j=0; !isnan(emp[i][j].val); j++)
            sum += emp[i][j].eff*emp[i][j].val;
       printf("i %d sum : %.2lE/%.2lE\n", i , sum, exp(2*gsl_sf_lnfact(i)+log(i+1)));
        for(j=0; !isnan(emp[i][j].val); j++) {
            emp[i][j].val = exp((emp[i][j].val)-2*gsl_sf_lnfact(i)-log((double)i+1.)+logNumberTrees(i+1));
            emp[i][j].eff /= numberTrees(i+1);
        }
        sum = 0;
        for(j=0; !isnan(emp[i][j].val); j++)
            sum += emp[i][j].eff*emp[i][j].val;
       printf("i %d sum : %lE\n", i , sum);
   }
}
