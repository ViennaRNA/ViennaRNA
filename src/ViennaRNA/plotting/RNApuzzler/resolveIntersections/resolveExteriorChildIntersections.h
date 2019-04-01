#ifndef RNAPUZZLER_RESOLVE_EXT_CHILD_INTERSECT
#define RNAPUZZLER_RESOLVE_EXT_CHILD_INTERSECT

#include "ViennaRNA/plotting/RNApuzzler/definitions.h"
#include "ViennaRNA/plotting/RNApuzzler/dataTypes/configtree_struct.h"
#include "ViennaRNA/plotting/RNApuzzler/dataTypes/tBaseInformation_struct.h"


PRIVATE void
resolveExteriorChildrenIntersectionXY(treeNode            *exteriorNode,
                                      short const *const  pair_table,
                                      const double        unpaired,
                                      const short         allowFlipping,
                                      double              *myX,
                                      double              *myY);


PRIVATE void
resolveExteriorChildrenIntersectionAffin(treeNode                 *exteriorNode,
                                         short const *const       pair_table,
                                         tBaseInformation *const  baseInformation,
                                         const double             unpaired,
                                         const short              allowFlipping);


PRIVATE void
resolveExteriorChildIntersections(treeNode                *exteriorNode,
                                  short const *const      pair_table,
                                  tBaseInformation *const baseInformation,
                                  const double            unpaired,
                                  const short             allowFlipping);


#include "resolveExteriorChildIntersections.inc"

#endif
