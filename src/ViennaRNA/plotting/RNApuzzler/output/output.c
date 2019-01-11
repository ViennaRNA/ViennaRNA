#include "ViennaRNA/RNApuzzler/output/output.h"
#include "ViennaRNA/RNApuzzler/vector_math.h"
#include "ViennaRNA/RNApuzzler/data/cfg_reader.h"
#include "ViennaRNA/RNApuzzler/data/configtree.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define COLOR_LIGHTGREEN "\x1b[92m"

#define COLORED_OUTPUT 0

/// Flags for special output
const short LOG_FLAG = 1;         /// print log messages

int printError(const char* fnName, const char * format, ... ) {
    FILE* stream = stdout;
    if (COLORED_OUTPUT) {
        fprintf(stream, ANSI_COLOR_RED);
    }
    va_list args;
    if (fnName != NULL) {
        fprintf(stream, "[ERROR] [%s] ", fnName);
    }
    va_start(args, format);
    int ret = vfprintf(stream, format, args);
    va_end(args);
    if (COLORED_OUTPUT) {
        fprintf(stream, ANSI_COLOR_RESET);
    }
    fflush(stream);
    return ret;
}

int printWarning(const char* fnName, const char * format, ... ) {
    FILE* stream = stdout;
    if (COLORED_OUTPUT) {
        fprintf(stream, ANSI_COLOR_YELLOW);
    }
    va_list args;
    if (fnName != NULL) {
        fprintf(stream, "[WARNING] [%s] ", fnName);
    }
    va_start(args, format);
    int ret = vfprintf(stream, format, args);
    va_end(args);
    if (COLORED_OUTPUT) {
        fprintf(stream, ANSI_COLOR_RESET);
    }
    fflush(stream);
    return ret;
}

int printInformation(const char* fnName, const char * format, ... ) {
    FILE* stream = stdout;
    if (COLORED_OUTPUT) {
        fprintf(stream, ANSI_COLOR_GREEN);
    }
    va_list args;
    if (fnName != NULL) {
        fprintf(stream, "[INFORMATION] [%s] ", fnName);
    }
    va_start(args, format);
    int ret = vfprintf(stream, format, args);
    va_end(args);
    if (COLORED_OUTPUT) {
        fprintf(stream, ANSI_COLOR_RESET);
    }
    fflush(stream);
    return ret;
}

int printDebug(const char* fnName, const char * format, ... ) {
    FILE* stream = stdout;
    if (COLORED_OUTPUT) {
        fprintf(stream, ANSI_COLOR_CYAN);
    }
    va_list args;
    if (fnName != NULL) {
        fprintf(stream, "[DEBUG] [%s] ", fnName);
    }
    va_start(args, format);
    int ret = vfprintf(stream, format, args);
    va_end(args);
    if (COLORED_OUTPUT) {
        fprintf(stream, ANSI_COLOR_RESET);
    }
    fflush(stream);
    return ret;
}

void printPath(
        char* fnName,
        treeNode** path,
        int pathLength,
        int nodeNumber
) {
    printDebug(fnName, "");
    printf("[PATH]");
    for (int i = 0; i < pathLength; i++) {
        if (COLORED_OUTPUT && i == nodeNumber) {
            printf(" "); printf(ANSI_COLOR_CYAN); printf("%c(%d)", getNodeName(path[i]), path[i]->childCount); printf(ANSI_COLOR_RESET);
        } else {
            printf(" %c(%d)", getNodeName(path[i]), path[i]->childCount);
        }
    }
    printf("\n");
}

void printConfigError(
    const char *fnName,
    const treeNode *node,
    const double *deltas
) {
    printError(fnName, "Invalid changes. How could that be?\n");

    config *cfg = node->cfg;
    int configSize = node->childCount + 1;

    for (int currentArc = 0; currentArc < configSize; ++currentArc) {
        double oldAngle = getArcAngle(cfg, currentArc);
        double newAngle = oldAngle + deltas[currentArc];
        short validAngle = 0 < newAngle && newAngle < MATH_TWO_PI;
        printError(fnName, "valid: %d %d[%d]: %+11.6lf° %+11.6lf° -> %+11.6lf°\n"
               , validAngle, getNodeID(node), currentArc
               , oldAngle, deltas[currentArc], newAngle);
    }
}

