#ifndef __VIENNA_RNA_PACKAGE_LFOLD_H__
#define __VIENNA_RNA_PACKAGE_LFOLD_H__

/**
 *  \file Lfold.h
 *  \brief Predicting local MFE structures of large sequences
 */

/**
 *  \brief The local analog to fold().
 * 
 *  Computes the minimum free energy structure including only base pairs
 *  with a span smaller than 'maxdist'
 * 
 *  \param string
 *  \param structure
 *  \param maxdist
 */
float Lfold(const char *string,
            char *structure,
            int maxdist);

/**
 *  \brief
 * 
 *  \param strings
 *  \param structure
 *  \param maxdist
 *  \return
 */
float aliLfold( const char **strings,
                char *structure,
                int maxdist);

/**
 *  \brief
 * 
 *  \param string
 *  \param structure
 *  \param maxdist
 *  \param zsc
 *  \param min_z
 */
float Lfoldz( const char *string,
              char *structure,
              int maxdist,
              int zsc,
              double min_z);

#endif
