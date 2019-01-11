#ifndef OUTPUT_H
#define OUTPUT_H

#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"

#include <stdarg.h>

void init_debug();
int printError(const char* fnName, const char * format, ... );
int printWarning(const char* fnName, const char * format, ... );
int printInformation(const char* fnName, const char * format, ... );
int printDebug(const char* fnName, const char * format, ... );

void printPath(
        char* fnName,
        treeNode** path,
        int pathLength,
        int nodeNumber
);

void printConfigError(
    const char *fnName,
    const treeNode *node,
    const double *deltas
);

#endif
