#ifndef VIENNA_RNA_PACKAGE_PROFILEALN_H
#define VIENNA_RNA_PACKAGE_PROFILEALN_H

float profile_aln(const float *T1,
                  const char *seq1,
                  const float *T2,
                  const char *seq2);

int set_paln_params(double gap_open,
                    double gap_ext,
                    double seqweight,
                    int free_ends);

#endif
