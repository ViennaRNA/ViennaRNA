#include "ViennaRNA/RNApuzzler/definitions.h"
#include "ViennaRNA/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/RNApuzzler/dataTypes/tBaseInformation_struct.h"


void resolveExteriorChildrenIntersectionXY(
        treeNode* exteriorNode,
        short const * const pair_table,
        const double unpaired,
        const short allowFlipping,
        double *myX,
        double *myY
);

void resolveExteriorChildrenIntersectionAffin(
        treeNode* exteriorNode,
        short const * const pair_table,
        tBaseInformation* const baseInformation,
        const double unpaired,
        const short allowFlipping
);

void resolveExteriorChildIntersections(
        treeNode* exteriorNode,
        short const * const pair_table,
        tBaseInformation* const baseInformation,
        const double unpaired,
        const short allowFlipping
);
