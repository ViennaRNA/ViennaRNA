#ifndef __VIENNA_RNA_PACKAGE_FIND_PATH_H__
#define __VIENNA_RNA_PACKAGE_FIND_PATH_H__

#include "data_structures.h"

/**
 *  \file findpath.h
 *  \brief Compute direct refolding paths between two secondary structures
 */

/**
 *  \brief Find energy of a saddle point between 2 structures
 *  (serch only direct path)
 *
 *  \param seq RNA sequence
 *  \param struc1 A pointer to the character array where the first
 *         secondary structure in dot-bracket notation will be written to
 *  \param struc2 A pointer to the character array where the second
 *         secondary structure in dot-bracket notation will be written to
 *  \param max integer how many strutures are being kept during the search
 *  \returns the saddle energy in 10cal/mol
 */
int     find_saddle(const char *seq,
                    const char *struc1,
                    const char *struc2,
                    int max);


/**
 *  \brief Find refolding path between 2 structures
 *  (serch only direct path)
 *
 *  \param seq RNA sequence
 *  \param s1 A pointer to the character array where the first
 *         secondary structure in dot-bracket notation will be written to
 *  \param s2 A pointer to the character array where the second
 *         secondary structure in dot-bracket notation will be written to
 *  \param maxkeep integer how many strutures are being kept during the search
 *  \returns direct refolding path between two structures
 */
path_t* get_path( const char *seq,
                  const char *s1,
                  const char* s2,
                  int maxkeep);

/**
 *  \brief Free memory allocated by get_path() function
 *
 *  \param path pointer to memory to be freed
 */
void    free_path(path_t *path);

#endif
