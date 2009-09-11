#ifndef __VIENNA_RNA_PACKAGE_STRING_DIST_H__
#define __VIENNA_RNA_PACKAGE_STRING_DIST_H__

#ifndef __VIENNA_RNA_PACKAGE_DIST_VARS_H__
#include "dist_vars.h"  /* defines the type Tree */
#endif
swString *Make_swString(char *string);
/* make input for string_edit_distance */
float     string_edit_distance(swString *T1, swString *T2);
/* compare two structures using string alignment */

#endif
