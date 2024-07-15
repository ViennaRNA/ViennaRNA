/** \file data_structures.c **/

/*
 *                Data structure creation/destruction
 *
 *                This file contains everything which is necessary to
 *                obtain and destroy datastructures used in the folding
 *                recurrences throughout the ViennaRNA package
 *
 *                c Ronny Lorenz
 *
 *                ViennaRNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ViennaRNA/utils/basic.h"

struct vrna_bt_stack_s {
  vrna_array(vrna_sect_t) stack;
};

/*
 #################################
 # PRIVATE MACROS                #
 #################################
 */

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

PUBLIC vrna_bts_t
vrna_bts_init(size_t  n)
{
  vrna_bts_t  bts = vrna_alloc(sizeof(struct vrna_bt_stack_s));

  if (n)
    vrna_array_init_size(bts->stack, n);
  else
    vrna_array_init(bts->stack)

  return bts;
}


PUBLIC void
vrna_bts_free(vrna_bts_t bts)
{
  if (bts) {
    vrna_array_free(bts->stack);
    free(bts);
  }
}


PUBLIC size_t
vrna_bts_push(vrna_bts_t  bts,
              vrna_sect_t element)
{
  if (bts) {
    vrna_array_append(bts->stack, element);
    return vrna_array_size(bts->stack);
  }
  
  return 0;
}


PUBLIC vrna_sect_t
vrna_bts_top(vrna_bts_t bts)
{
  if ((bts) &&
      (vrna_array_size(bts->stack) > 0)) {
    return bts->stack[vrna_array_size(bts->stack) - 1];
  }
  
  return (vrna_sect_t){0};
}


PUBLIC vrna_sect_t
vrna_bts_pop(vrna_bts_t bts)
{
  if ((bts) &&
      (vrna_array_size(bts->stack) > 0)) {
    return bts->stack[--(vrna_array_size(bts->stack))];
  }
  
  return (vrna_sect_t){0};
}


PUBLIC size_t
vrna_bts_size(vrna_bts_t bts)
{
  if (bts)
    return vrna_array_size(bts->stack);
  
  return 0;
}


#ifndef VRNA_DISABLE_C11_FEATURES
PUBLIC void
vrna_C11_features(void){}
#endif
