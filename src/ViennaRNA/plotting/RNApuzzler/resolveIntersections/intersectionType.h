#ifndef INTERSECTION_TYPE_H
#define INTERSECTION_TYPE_H

typedef enum {
  noIntersection = 0,
  LxL = 1,
  LxS = 2,
  SxL = 3,
  SxS = 4,
  LxB = 5,
  BxL = 6,
  SxB = 7,
  BxS = 8,
  BxB = 9,
  siblings = 10,
  exterior = 11,
} intersectionType;

char *intersectionTypeToString(
    const intersectionType it
);

#endif
