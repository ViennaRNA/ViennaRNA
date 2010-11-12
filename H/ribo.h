#ifndef __VIENNA_RNA_PACKAGE_RIBOSUM_H__
#define __VIENNA_RNA_PACKAGE_RIBOSUM_H__

float **get_ribosum(const char **Alseq, int n_seq, int length);
float **get_ribosum_slice(const char **Alseq, int n_seq, int start, int length);

#endif
