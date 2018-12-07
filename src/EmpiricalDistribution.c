#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "Utils.h"
#include "EmpiricalDistribution.h"



int compareEmpiricalDistrivutionItem(const void* a, const void* b) {
    if(((TypeEmpiricalDistributionItem*)a)->val>((TypeEmpiricalDistributionItem*)b)->val)
        return 1;
    if(((TypeEmpiricalDistributionItem*)a)->val<((TypeEmpiricalDistributionItem*)b)->val)
        return -1;
    return 0;
}

