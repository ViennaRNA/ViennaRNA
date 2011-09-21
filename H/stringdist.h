#ifndef __VIENNA_RNA_PACKAGE_STRING_DIST_H__
#define __VIENNA_RNA_PACKAGE_STRING_DIST_H__

/**
 *  \file stringdist.h
 *  \brief Functions for String Alignment
 */

#include "dist_vars.h"


/**
 *  \brief Convert a structure into a format suitable for string_edit_distance().
 * 
 *  \param string
 *  \return
 */
swString *Make_swString(char *string);

/**
 *  \brief Calculate the string edit distance of T1 and T2.
 * 
 *  \param  T1
 *  \param  T2
 *  \return
 */
float     string_edit_distance( swString *T1,
                                swString *T2);

#endif
