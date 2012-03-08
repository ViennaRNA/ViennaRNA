#ifndef __VIENNA_RNA_PACKAGE_GQUAD_H__
#define __VIENNA_RNA_PACKAGE_GQUAD_H__

#include "data_structures.h"

/**
 *  \file gquad.h
 *  \brief Various functions related to G-quadruplex computations
 */

#define   VRNA_GQUAD_MAX_STACK_SIZE     7
#define   VRNA_GQUAD_MIN_STACK_SIZE     2
#define   VRNA_GQUAD_MAX_LINKER_LENGTH  25
#define   VRNA_GQUAD_MIN_LINKER_LENGTH  1
#define   VRNA_GQUAD_MIN_BOX_SIZE       (4*VRNA_GQUAD_MIN_STACK_SIZE)+(3*VRNA_GQUAD_MIN_LINKER_LENGTH)
#define   VRNA_GQUAD_MAX_BOX_SIZE       (4*VRNA_GQUAD_MAX_STACK_SIZE)+(3*VRNA_GQUAD_MAX_LINKER_LENGTH)
#define   VRNA_GQUAD_MISMATCH_PENALTY   300   /* penalty for incompatible nucleotides in an alignment that destruct a gquad layer */
#define   VRNA_GQUAD_MISMATCH_NUM_ALI   1   /* maximum number of mismatching sequences in the alignment when gquad should be formed */


int         gquad_contribution( int L,
                                int l1,
                                int l2,
                                int l3);

int         gquad_ali_contribution( int i,
                                    int L,
                                    int l1,
                                    int l2,
                                    int l3,
                                    short **S,
                                    int n_seq);

void        gquad_ali_contribution_en(int i,
                                      int L,
                                      int l1,
                                      int l2,
                                      int l3,
                                      const short **S,
                                      int n_seq,
                                      int en[2]);

int         *get_gquad_matrix(short *S);

int         *get_gquad_ali_matrix(short *S_cons,
                                  short **S,
                                  int n_seq,
                                  int *pscore);

FLT_OR_DBL  *get_gquad_pf_matrix( short *S,
                                  FLT_OR_DBL *scale);

plist       *get_plist_gquad_from_pr( short *S,
                                      int gi,
                                      int gj,
                                      FLT_OR_DBL *G,
                                      FLT_OR_DBL *probs,
                                      FLT_OR_DBL *scale);
plist       *get_plist_gquad_from_pr_max(short *S,
                                      int gi,
                                      int gj,
                                      FLT_OR_DBL *G,
                                      FLT_OR_DBL *probs,
                                      FLT_OR_DBL *scale,
                                      int *L,
                                      int l[3]);

plist       *get_plist_gquad_from_db( const char *structure,
                                      float pr);

plist       *Gquadcomputeinnerprobability(short *S,
                                          FLT_OR_DBL *G,
                                          FLT_OR_DBL  *probs,
                                          FLT_OR_DBL *scale);

int         parse_gquad(const char *struc, int *L, int l[3]);



#endif
