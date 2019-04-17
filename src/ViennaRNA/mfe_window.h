#ifndef VIENNA_RNA_PACKAGE_MFE_WINDOW_H
#define VIENNA_RNA_PACKAGE_MFE_WINDOW_H

#include <stdio.h>
#include <ViennaRNA/fold_compound.h>

/**
 *
 *  @file mfe_window.h
 *  @ingroup  mfe, mfe_window
 *  @brief Compute local Minimum Free Energy (MFE) using a sliding window approach and
 *         backtrace corresponding secondary structures
 *
 *  This file includes the interface to all functions related to predicting locally
 *  stable secondary structures.
 */


/**
 *  @addtogroup  mfe_window
 *  @{
 *  @brief  Variations of the local (sliding window) Minimum Free Energy (MFE) prediction algorithm
 *
 *  We provide implementations for the local (sliding window) MFE prediction algorithm for
 *  * Single sequences,
 *  * Multiple sequence alignments (MSA), and
 *
 *  Note, that our implementation scans an RNA sequence (or MSA) from the 3' to the 5'
 *  end, and reports back locally optimal (consensus) structures, the corresponding free
 *  energy, and the position of the sliding window in global coordinates.
 *
 *  For any particular RNA sequence (or MSA) multiple locally optimal (consensus)
 *  secondary structures may be predicted. Thus, we tried to implement an interface that
 *  allows for an effortless conversion of the corresponding hits into any target data
 *  structure. As a consequence, we provide two distinct ways to retrieve the corresponding
 *  predictions, either
 *  * through directly writing to an open @p FILE stream on-the-fly, or
 *  * through a callback function mechanism.
 *
 *  The latter allows one to store the results in any possible target data structure. Our
 *  implementations then pass the results through the user-implemented callback as soon as
 *  the prediction for a particular window is finished.
 */

/**
 *  @brief  The default callback for sliding window MFE structure predictions
 *
 *  @callback
 *  @parblock
 *  This function will be called for each hit in a sliding window MFE prediction.
 *  @endparblock
 *  @see vrna_mfe_window()
 *
 *  @param start provides the first position of the hit (1-based, relative to entire sequence/alignment)
 *  @param end provides the last position of the hit (1-based, relative to the entire sequence/alignment)
 *  @param structure provides the (sub)structure in dot-bracket notation
 *  @param en is the free energy of the structure hit in kcal/mol
 *  @param data is some arbitrary data pointer passed through by the function executing the callback
 */
typedef void (vrna_mfe_window_callback)(int         start,
                                        int         end,
                                        const char  *structure,
                                        float       en,
                                        void        *data);


#ifdef VRNA_WITH_SVM
typedef void (vrna_mfe_window_zscore_callback)(int        start,
                                               int        end,
                                               const char *structure,
                                               float      en,
                                               float      zscore,
                                               void       *data);
#endif

/**
 *  @name Basic local (sliding window) MFE prediction interface
 *  @{
 */

/**
 *  @brief Local MFE prediction using a sliding window approach.
 *
 *  Computes minimum free energy structures using a sliding window
 *  approach, where base pairs may not span outside the window.
 *  In contrast to vrna_mfe(), where a maximum base pair span
 *  may be set using the #vrna_md_t.max_bp_span attribute and one
 *  globally optimal structure is predicted, this function uses a
 *  sliding window to retrieve all locally optimal structures within
 *  each window.
 *  The size of the sliding window is set in the #vrna_md_t.window_size
 *  attribute, prior to the retrieval of the #vrna_fold_compound_t
 *  using vrna_fold_compound() with option #VRNA_OPTION_WINDOW
 *
 *  The predicted structures are written on-the-fly, either to
 *  stdout, if a NULL pointer is passed as file parameter, or to
 *  the corresponding filehandle.
 *
 *  @see  vrna_fold_compound(), vrna_mfe_window_zscore(), vrna_mfe(),
 *        vrna_Lfold(), vrna_Lfoldz(),
 *        #VRNA_OPTION_WINDOW, #vrna_md_t.max_bp_span, #vrna_md_t.window_size
 *
 *  @param  vc        The #vrna_fold_compound_t with preallocated memory for the DP matrices
 *  @param  file      The output file handle where predictions are written to (maybe NULL)
 */
float
vrna_mfe_window(vrna_fold_compound_t  *vc,
                FILE                  *file);


float
vrna_mfe_window_cb(vrna_fold_compound_t     *vc,
                   vrna_mfe_window_callback *cb,
                   void                     *data);


#ifdef VRNA_WITH_SVM
/**
 *  @brief Local MFE prediction using a sliding window approach (with z-score cut-off)
 *
 *  Computes minimum free energy structures using a sliding window
 *  approach, where base pairs may not span outside the window.
 *  This function is the z-score version of vrna_mfe_window(), i.e.
 *  only predictions above a certain z-score cut-off value are
 *  printed.
 *  As for vrna_mfe_window(), the size of the sliding window is set in
 *  the #vrna_md_t.window_size attribute, prior to the retrieval of
 *  the #vrna_fold_compound_t using vrna_fold_compound() with option
 *  #VRNA_OPTION_WINDOW.
 *
 *  The predicted structures are written on-the-fly, either to
 *  stdout, if a NULL pointer is passed as file parameter, or to
 *  the corresponding filehandle.
 *
 *  @see  vrna_fold_compound(), vrna_mfe_window_zscore(), vrna_mfe(),
 *        vrna_Lfold(), vrna_Lfoldz(),
 *        #VRNA_OPTION_WINDOW, #vrna_md_t.max_bp_span, #vrna_md_t.window_size
 *
 *  @param  vc        The #vrna_fold_compound_t with preallocated memory for the DP matrices
 *  @param  min_z     The minimal z-score for a predicted structure to appear in the output
 *  @param  file      The output file handle where predictions are written to (maybe NULL)
 */
float
vrna_mfe_window_zscore(vrna_fold_compound_t *vc,
                       double               min_z,
                       FILE                 *file);


float
vrna_mfe_window_zscore_cb(vrna_fold_compound_t            *vc,
                          double                          min_z,
                          vrna_mfe_window_zscore_callback *cb,
                          void                            *data);


#endif

/* End basic local MFE interface */
/**@}*/

/**
 *  @name Simplified local MFE prediction using sequence(s) or multiple sequence alignment(s)
 *  @{
 */

/**
 *  @brief Local MFE prediction using a sliding window approach (simplified interface)
 *
 *  This simplified interface to vrna_mfe_window() computes the MFE and locally
 *  optimal secondary structure using default options. Structures are predicted
 *  using a sliding window approach, where base pairs may not span outside the
 *  window. Memory required for dynamic programming (DP) matrices will
 *  be allocated and free'd on-the-fly. Hence, after return of this function, the recursively filled
 *  matrices are not available any more for any post-processing.
 *
 *  @note In case you want to use the filled DP matrices for any subsequent post-processing step, or
 *  you require other conditions than specified by the default model details, use vrna_mfe_window(),
 *  and the data structure #vrna_fold_compound_t instead.
 *
 *  @see  vrna_mfe_window(), vrna_Lfoldz(), vrna_mfe_window_zscore()
 *
 *  @param  string      The nucleic acid sequence
 *  @param  window_size The window size for locally optimal structures
 *  @param  file        The output file handle where predictions are written to (if NULL, output is written to stdout)
 */
float
vrna_Lfold(const char *string,
           int        window_size,
           FILE       *file);


float
vrna_Lfold_cb(const char                *string,
              int                       window_size,
              vrna_mfe_window_callback  *cb,
              void                      *data);


#ifdef VRNA_WITH_SVM
/**
 *  @brief Local MFE prediction using a sliding window approach with z-score cut-off (simplified interface)
 *
 *  This simplified interface to vrna_mfe_window_zscore() computes the MFE and locally
 *  optimal secondary structure using default options. Structures are predicted
 *  using a sliding window approach, where base pairs may not span outside the
 *  window. Memory required for dynamic programming (DP) matrices will
 *  be allocated and free'd on-the-fly. Hence, after return of this function, the recursively filled
 *  matrices are not available any more for any post-processing.
 *  This function is the z-score version of vrna_Lfold(), i.e.
 *  only predictions above a certain z-score cut-off value are
 *  printed.
 *
 *  @note In case you want to use the filled DP matrices for any subsequent post-processing step, or
 *  you require other conditions than specified by the default model details, use vrna_mfe_window(),
 *  and the data structure #vrna_fold_compound_t instead.
 *
 *  @see  vrna_mfe_window_zscore(), vrna_Lfold(), vrna_mfe_window()
 *
 *  @param  string      The nucleic acid sequence
 *  @param  window_size The window size for locally optimal structures
 *  @param  min_z       The minimal z-score for a predicted structure to appear in the output
 *  @param  file        The output file handle where predictions are written to (if NULL, output is written to stdout)
 */
float
vrna_Lfoldz(const char  *string,
            int         window_size,
            double      min_z,
            FILE        *file);


float
vrna_Lfoldz_cb(const char                       *string,
               int                              window_size,
               double                           min_z,
               vrna_mfe_window_zscore_callback  *cb,
               void                             *data);


#endif

float vrna_aliLfold(const char  **alignment,
                    int         maxdist,
                    FILE        *fp);


float vrna_aliLfold_cb(const char               **alignment,
                       int                      maxdist,
                       vrna_mfe_window_callback *cb,
                       void                     *data);


/* End simplified local MFE interface */
/**@}*/

/* End group mfe_fold_window */
/**@}*/


#endif
