#include "ViennaRNA/RNApuzzler/resolveIntersections/intersectionType.h"

char *intersectionTypeToString(
    const intersectionType it
) {
      switch (it) {
        case LxL:      return "LxL";
        case LxS:      return "LxS";
        case LxB:      return "LxB";
        case SxL:      return "SxL";
        case SxS:      return "SxS";
        case SxB:      return "SxB";
        case BxL:      return "BxL";
        case BxS:      return "BxS";
        case BxB:      return "BxB";
        case siblings: return "BRA";
        case exterior: return "EXT";
        default:       return "UNK";
    }
}


