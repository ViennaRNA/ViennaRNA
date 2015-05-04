#ifndef VIENNA_RNA_PACKAGE_TREE_DIST_H
#define VIENNA_RNA_PACKAGE_TREE_DIST_H

/**
 *  \file treedist.h
 *  \brief Functions for Tree Edit Distances
 */

#include <ViennaRNA/dist_vars.h>

/**
 *  \brief Constructs a Tree ( essentially the postorder list ) of the
 *  structure 'struc', for use in tree_edit_distance().
 * 
 *  \param  struc may be any rooted structure representation.
 *  \return
 */
Tree   *make_tree(char *struc);

/**
 *  \brief Calculates the edit distance of the two trees.
 * 
 *  \param T1
 *  \param T2
 *  \return
 */
float   tree_edit_distance( Tree *T1,
                            Tree *T2);

/**
 *  \brief Print a tree (mainly for debugging)
 */
void    print_tree(Tree *t);

/**
 *  \brief Free the memory allocated for Tree t.
 * 
 *  \param t
 */
void    free_tree(Tree *t);

#endif
