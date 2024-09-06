/*             Global variables for Distance-Package */
#include "ViennaRNA/dist_vars.h"

int   edit_backtrack = 0; /* calculate aligned representation */

char  *aligned_line[4];   /* containes the aligned string representations */

int   cost_matrix = 0;    /* 0 for usual costs, 1 for Shapiro's costs */
