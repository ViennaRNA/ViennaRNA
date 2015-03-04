#ifndef VIENNA_RNA_PACKAGE_RIBOSUM_H
#define VIENNA_RNA_PACKAGE_RIBOSUM_H

float **get_ribosum(const char **Alseq,
                    int n_seq,
                    int length);

/**
 *  \brief Read a ribosum or other user-defined scoring matrix
 * 
 *  \ingroup consensus_fold
 * 
 */
float   **readribosum(char *name);

#endif
