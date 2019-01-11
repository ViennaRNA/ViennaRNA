#ifndef GEOGEBRA_OUTPUT_H
#define GEOGEBRA_OUTPUT_H

#include "ViennaRNA/RNApuzzler/dataTypes/boundingBoxes_struct.h"
#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/resolveIntersections/intersectionType.h"

void GEOGEBRA_printStem(stemBox* rect);
void GEOGEBRA_printLoop(loopBox* circ);
void GEOGEBRA_printNode(treeNode* node, short printStem, short printLoop);
void GEOGEBRA_printPath(treeNode* nodeTOP, treeNode* nodeBOTTOM);
void GEOGEBRA_printTree(treeNode* node);
void GEOGEBRA_generatePath(treeNode* nodeTOP, treeNode* nodeBOTTOM, int changeNumber);
void GEOGEBRA_generateTree(treeNode* node, int changeNumber);
void GEOGEBRA_generateTreeWithIntersectionType(treeNode* node, int changeNumber, intersectionType it, short step);

void GEOGEBRA_printLxL(treeNode* source, treeNode* target);
void GEOGEBRA_printLxS(treeNode* source, treeNode* target);
void GEOGEBRA_printSxL(treeNode* source, treeNode* target);
void GEOGEBRA_printSxS(treeNode* source, treeNode* target);

void GEOGEBRA_printCircle(const char* name, const double c[2], const double r);
void GEOGEBRA_printLinePointPoint(const char* name, const double p1[2], const double p2[2]);
void GEOGEBRA_printLinePointDir(const char* name, const double p[2], const double v[2]);
void GEOGEBRA_printRay(const char* name, const double p[2], const double v[2]);
void GEOGEBRA_printPoint(const char* name, const double p[2]);
void GEOGEBRA_printRectangle(const char* name, const double a[2], const double b[2], const double c[2], const double e[2]);

void GEOGEBRA_printInteractiveTree(const treeNode* tree);

#endif
