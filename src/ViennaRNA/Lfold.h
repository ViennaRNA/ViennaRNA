#ifndef VIENNA_RNA_PACKAGE_LFOLD_H
#define VIENNA_RNA_PACKAGE_LFOLD_H

/**
 */

/* make this interface backward compatible with RNAlib < 2.2.0 */
#define VRNA_BACKWARD_COMPAT

#ifdef __GNUC__
#define DEPRECATED(func) func __attribute__ ((deprecated))
#else
#define DEPRECATED(func) func
#endif

/**
 *  @addtogroup local_fold
 *
 *  Local structures can be predicted by a modified version of the
 *  fold() algorithm that restricts the span of all base pairs.
 *  @{
 *    @file Lfold.h
 *    @brief Predicting local MFE structures of large sequences
 *
 *  @}
 */

/**
 *  @addtogroup local_mfe_fold
 *  @{
 *
 *  @}
 */

/**
 *  @brief Local MFE prediction using a sliding window approach.
 *
 *  Computes minimum free energy structures using a sliding window
 *  approach, where base pairs may not span outside the window.
 *  In contrast to vrna_fold(), where a maximum base pair span
 *  may be set using the #vrna_md_t.max_bp_span attribute and one
 *  globally optimal structure is predicted, this function uses a
 *  sliding window to retrieve all locally optimal structures within
 *  each window.
 *  The size of the sliding window is set in the #vrna_md_t.window_size
 *  attribute, prior to the retrieval of the #vrna_fold_compound
 *  using vrna_get_fold_compound() with option #VRNA_OPTION_WINDOW
 *
 *  The predicted structures are written on-the-fly, either to
 *  stdout, if a NULL pointer is passed as file parameter, or to
 *  the corresponding filehandle.
 *
 *  @ingroup local_mfe_fold
 * 
 *  @see  vrna_get_fold_compound(), vrna_Lfoldz(), vrna_fold(),
 *        #VRNA_OPTION_WINDOW, #vrna_md_t.max_bp_span, #vrna_md_t.window_size
 *
 *  @param  vc        The #vrna_fold_compound with preallocated memory for the DP matrices
 *  @param  file      The output file handle where predictions are written to (maybe NULL)
 */
float vrna_Lfold( vrna_fold_compound *vc, FILE *file);

#ifdef USE_SVM
/**
 *  @brief Local MFE prediction using a sliding window approach (with z-score cut-off)
 *
 *  Computes minimum free energy structures using a sliding window
 *  approach, where base pairs may not span outside the window.
 *  This function is the z-score version of vrna_Lfold(), i.e.
 *  only predictions above a certain z-score cut-off value are
 *  printed.
 *  As for vrna_Lfold(), the size of the sliding window is set in
 *  the #vrna_md_t.window_size attribute, prior to the retrieval of
 *  the #vrna_fold_compound using vrna_get_fold_compound() with option
 *  #VRNA_OPTION_WINDOW.
 *
 *  The predicted structures are written on-the-fly, either to
 *  stdout, if a NULL pointer is passed as file parameter, or to
 *  the corresponding filehandle.
 *
 *  @ingroup local_mfe_fold
 * 
 *  @see  vrna_get_fold_compound(), vrna_Lfoldz(), vrna_fold(),
 *        #VRNA_OPTION_WINDOW, #vrna_md_t.max_bp_span, #vrna_md_t.window_size
 *
 *  @param  vc        The #vrna_fold_compound with preallocated memory for the DP matrices
 *  @param  min_z     The minimal z-score for a predicted structure to appear in the output
 *  @param  file      The output file handle where predictions are written to (maybe NULL)
 */
float vrna_Lfoldz(vrna_fold_compound *vc, double min_z, FILE *file);
#endif


/**
 *  @addtogroup local_consensus_fold
 *  @{
 *
 *  @}
 */

/**
 *  @brief
 *
 *  @ingroup local_consensus_fold
 * 
 *  @param strings
 *  @param structure
 *  @param maxdist
 *  @return
 */
float aliLfold( const char **strings,
                char *structure,
                int maxdist);

#ifdef  VRNA_BACKWARD_COMPAT

/**
 *  @brief The local analog to fold().
 * 
 *  Computes the minimum free energy structure including only base pairs
 *  with a span smaller than 'maxdist'
 *
 *  @ingroup local_mfe_fold
 *
 *  @deprecated Use vrna_Lfold() instead!
 */
DEPRECATED(float Lfold(const char *string, char *structure, int maxdist));

/**
 *  @brief
 * 
 *  @ingroup local_mfe_fold
 * 
 *  @deprecated Use vrna_Lfoldz() instead!
 */
DEPRECATED(float Lfoldz(const char *string, char *structure, int maxdist, int zsc, double min_z));

#endif

#endif
