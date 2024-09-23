#ifndef VIENNA_RNA_PACKAGE_STRUCTURES_HELIX_H
#define VIENNA_RNA_PACKAGE_STRUCTURES_HELIX_H

/**
 *  @file     ViennaRNA/structures/helix.h
 *  @ingroup  struct_utils
 *  @brief    Functions for secondary structure helices
 */


/**
 *  @addtogroup struct_utils_helix_list
 *  @{
 */

/**
 *  @brief Convenience typedef for data structure #vrna_hx_s
 *  @ingroup  struct_utils_helix_list
 */
typedef struct vrna_hx_s vrna_hx_t;


/**
 *  @brief  Data structure representing an entry of a helix list
 */
struct vrna_hx_s {
  unsigned int  start;
  unsigned int  end;
  unsigned int  length;
  unsigned int  up5;
  unsigned int  up3;
};


/**
 *  @brief  Convert a pair table representation of a secondary structure into a helix list
 *
 *  @param  pt  The secondary structure in pair table representation
 *  @return     The secondary structure represented as a helix list
 */
vrna_hx_t *
vrna_hx_from_ptable(short *pt);


/**
 *  @brief  Create a merged helix list from another helix list
 */
vrna_hx_t *
vrna_hx_merge(const vrna_hx_t *list,
              int             maxdist);


/* End helix list interface */
/** @} */

#endif
