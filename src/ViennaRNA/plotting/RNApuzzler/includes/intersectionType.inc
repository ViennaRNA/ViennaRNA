#ifndef RNAPUZZLER_INTERSECTION_TYPE_H
#define RNAPUZZLER_INTERSECTION_TYPE_H

/*
 *      RNApuzzler intersection type enum
 *
 *      c  Daniel Wiegreffe, Daniel Alexander, Dirk Zeckzer
 *      ViennaRNA package
 */

typedef enum {
  noIntersection  = 0,
  LxL             = 1,
  LxS             = 2,
  SxL             = 3,
  SxS             = 4,
  LxB             = 5,
  BxL             = 6,
  SxB             = 7,
  BxS             = 8,
  BxB             = 9,
  siblings        = 10,
  exterior        = 11,
} intersectionType;


PRIVATE char *
intersectionTypeToString(const intersectionType it)
{
  switch (it) {
    case LxL:
      return "LxL";
    case LxS:
      return "LxS";
    case LxB:
      return "LxB";
    case SxL:
      return "SxL";
    case SxS:
      return "SxS";
    case SxB:
      return "SxB";
    case BxL:
      return "BxL";
    case BxS:
      return "BxS";
    case BxB:
      return "BxB";
    case siblings:
      return "BRA";
    case exterior:
      return "EXT";
    default:
      return "UNK";
  }
}


#endif
