#ifndef VIENNA_RNA_PACKAGE_RIBOSUM_H
#define VIENNA_RNA_PACKAGE_RIBOSUM_H

/**
 *  @addtogroup   file_utils
 *
 *  @{
 *
 *  @file ribo.h
 *  @brief  Parse RiboSum Scoring Matrices for Covariance Scoring of Alignments
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
