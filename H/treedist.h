#ifndef DIST_VARS_H
#include "dist_vars.h"  /* defines the type Tree */
#endif
extern  Tree   *make_tree(char *struc); /* make input for tree_edit_distance */
extern  float   tree_edit_distance(Tree *T1, Tree *T2);
/* compare to structures using tree editing */
extern  void    print_tree(Tree *t);    /* mainly for debugging */
extern  void    free_tree(Tree *t);     /* free space allocated by make_tree */
