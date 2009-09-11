#ifndef __VIENNA_RNA_PACKAGE_TREE_DIST_H__
#define __VIENNA_RNA_PACKAGE_TREE_DIST_H__

#ifndef __VIENNA_RNA_PACKAGE_DIST_VARS_H__
#include "dist_vars.h"  /* defines the type Tree */
#endif
Tree   *make_tree(char *struc); /* make input for tree_edit_distance */
float   tree_edit_distance(Tree *T1, Tree *T2);
/* compare two structures using tree editing */
void    print_tree(Tree *t);    /* mainly for debugging */
void    free_tree(Tree *t);     /* free space allocated by make_tree */
#endif
