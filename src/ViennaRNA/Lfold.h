#ifndef VIENNA_RNA_PACKAGE_LFOLD_H
#define VIENNA_RNA_PACKAGE_LFOLD_H

/**
 */

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

/**
 *  \addtogroup local_fold
 *
 *  Local structures can be predicted by a modified version of the
 *  fold() algorithm that restricts the span of all base pairs.
 *  @{
 *    \file Lfold.h
 *    \brief Predicting local MFE structures of large sequences
 *
 *  @}
 */

/**
 *  \addtogroup local_mfe_fold
 *  @{
 *
 *  @}
 */

/**
 *  \brief Local MFE prediction using a sliding window approach.
 *
 *  Computes minimum free energy structures using a sliding window
 *  approach, where base pairs may not span outside the window.
 *  In contrast to vrna_fold(), where the maximum base pair span
 *  may be set using #vrna_md_t.max_bp_span and one globally optimal
 *  structure is predicted, this function uses a sliding window to
 *  retrieve all locally optimal structures within each window.
 *  The predicted structures are written on-the-fly, either to
 *  stdout, if a NULL pointer is passed as file parameter, or to
 *  the corresponding filehandle.
 *
 *  \ingroup local_mfe_fold
 * 
 *  @see vrna_get_fold_compound_window(), vrna_Lfoldz(), vrna_fold
 *
 *  \param vc   
 *  \param structure
 *  \param maxdist
 */
float vrna_Lfold( vrna_fold_compound *vc, FILE *file);

#ifdef USE_SVM
float vrna_Lfoldz(vrna_fold_compound *vc, double min_z, FILE *file);
#endif


/**
 *  \addtogroup local_consensus_fold
 *  @{
 *
 *  @}
 */

/**
 *  \brief
 *
 *  \ingroup local_consensus_fold
 * 
 *  \param strings
 *  \param structure
 *  \param maxdist
 *  \return
 */
float aliLfold( const char **strings,
                char *structure,
                int maxdist);

#ifdef  VRNA_BACKWARD_COMPAT

/**
 *  \brief The local analog to fold().
 * 
 *  Computes the minimum free energy structure including only base pairs
 *  with a span smaller than 'maxdist'
 *
 *  \ingroup local_mfe_fold
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
 *  \ingroup local_mfe_fold
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

#endif
