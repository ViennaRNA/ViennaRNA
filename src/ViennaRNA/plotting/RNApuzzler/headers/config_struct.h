#ifndef RNAPUZZLER_CONFIG_STRUCT_H
#define RNAPUZZLER_CONFIG_STRUCT_H

#include "configArc_struct.h"

/**
 * @brief The config struct implements a linked list to store configurations for a single loop.
 * Iteration through all arcs is possible using the next-member.
 */
typedef struct {
  // current radius
  double    radius;

  // optimal radius for current status
  double    minRadius;

  // default radius
  double    defaultRadius;

  // arcs
  configArc *cfgArcs;
  int       numberOfArcs;
} config;

#endif
