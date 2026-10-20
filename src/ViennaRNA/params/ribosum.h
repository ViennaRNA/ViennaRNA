#ifndef VIENNA_RNA_PACKAGE_RIBOSUM_H
#define VIENNA_RNA_PACKAGE_RIBOSUM_H

/**
 *  @file ViennaRNA/params/ribosum.h
 *  @ingroup   file_utils
 *  @brief  Parse RiboSum Scoring Matrices for Covariance Scoring of Alignments
 */


double **
vrna_ribosum_mx(double min_identity_pct,
                double max_identity_pct);


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup   file_utils
 *  @{
 */

/**
 *  @brief Retrieve a RiboSum Scoring Matrix for a given Alignment
 *  @ingroup energy_parameters
 *
 */
float **get_ribosum(const char  **Alseq,
                    int         n_seq,
                    int         length);


/**
 *  @brief Read a RiboSum or other user-defined Scoring Matrix and Store into global Memory
 *
 */
float **readribosum(char *name);


/**
 *  @}
 */
#endif

#endif
