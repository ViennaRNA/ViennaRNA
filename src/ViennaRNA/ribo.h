#ifndef VIENNA_RNA_PACKAGE_RIBOSUM_H
#define VIENNA_RNA_PACKAGE_RIBOSUM_H

/**
 *  @file ribo.h
 *  @ingroup   file_utils
 *  @brief  Parse RiboSum Scoring Matrices for Covariance Scoring of Alignments
 */

/**
 *  @{
 *  @ingroup   file_utils
 */

/**
 *  @brief Retrieve a RiboSum Scoring Matrix for a given Alignment
 *  \ingroup consensus_fold
 * 
 */
float **get_ribosum(const char **Alseq,
                    int n_seq,
                    int length);

/**
 *  \brief Read a RiboSum or other user-defined Scoring Matrix and Store into global Memory
 * 
 */
float   **readribosum(char *name);

/**
 *  @}
 */
#endif
