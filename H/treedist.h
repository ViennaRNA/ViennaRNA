#ifndef DIST_VARS_H
#include "dist_vars.h"  /* defines the type Tree */
#endif
extern  Tree   *make_tree(char *struc);
extern  float   tree_edit_distance(Tree *T1, Tree *T2);
extern  void    print_tree(Tree *t);
extern  void    free_tree(Tree *t);
