#ifndef RNAPUZZLER_STRUCTURE_STRUCT_H
#define RNAPUZZLER_STRUCTURE_STRUCT_H

#include "config_struct.h"

typedef enum {
  TYPE_BASE_NONE  = 0,
  TYPE_EXTERIOR   = 1,
  TYPE_STEM       = 2,
  TYPE_BULGE      = 3,
  TYPE_LOOP1      = 4,
  TYPE_LOOP2      = 5
} BASE_TYPE;

typedef struct {
  BASE_TYPE baseType;
  double    angle;
  double    distance;
  config    *config;
} tBaseInformation;

#endif
