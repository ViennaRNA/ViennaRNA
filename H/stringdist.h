#ifndef DIST_VARS_H
#include "dist_vars.h"  /* defines the type Tree */
#endif
extern swString *Make_swString(char *string);
/* make input for string_edit_distance */
extern float     string_edit_distance(swString *T1, swString *T2);
/* compare to structures using string alignment */
