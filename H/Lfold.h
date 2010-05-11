#ifndef __VIENNA_RNA_PACKAGE_LFOLD_H__
#define __VIENNA_RNA_PACKAGE_LFOLD_H__

float  Lfold(const char *string, char *structure, int maxdist);
float  aliLfold(char **strings, char *structure, int maxdist);

#endif
