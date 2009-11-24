#ifndef __VIENNA_RNA_PACKAGE_FIND_PATH_H__
#define __VIENNA_RNA_PACKAGE_FIND_PATH_H__

#include "data_structures.h"

int     find_saddle (char *seq, char *struc1, char *struc2, int max);
path_t* get_path(char *seq, char *s1, char* s2, int maxkeep);

#endif
