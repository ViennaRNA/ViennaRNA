#ifndef __VIENNA_RNA_PACKAGE_GQUAD_H__
#define __VIENNA_RNA_PACKAGE_GQUAD_H__

#include "data_structures.h"

/**
 *  \file gquad.h
 *  \brief Various functions related to G-quadruplex computations
 */


int         E_gquad(int L,
                    int l[3],
                    paramT *P);

FLT_OR_DBL exp_E_gquad( int L,
                        int l[3],
                        pf_paramT *pf);

int         E_gquad_ali(int i,
                        int L,
                        int l[3],
                        const short **S,
                        int n_seq,
                        paramT *P);


void        E_gquad_ali_en( int i,
                            int L,
                            int l[3],
                            const short **S,
                            int n_seq,
                            int en[2],
                            paramT *P);

/**
 *  \brief Get a triangular matrix prefilled with minimum free energy
 *  contributions of G-quadruplexes.
 *
 *  At each position ij in the matrix, the minimum free energy of any
 *  G-quadruplex delimited by i and j is stored. If no G-quadruplex formation
 *  is possible, the matrix element is set to INF.
 *  Access the elements in the matrix via matrix[indx[j]+i]. To get
 *  the integer array indx see get_jindx().
 *
 *  \see get_jindx(), encode_sequence()
 *
 *  \param S  The encoded sequence
 *  \param P  A pointer to the data structure containing the precomputed energy contributions
 *  \return   A pointer to the G-quadruplex contribution matrix
*/
int         *get_gquad_matrix(short *S, paramT *P);

int         *get_gquad_ali_matrix(short *S_cons,
                                  short **S,
                                  int n_seq,
                                  paramT *P);

FLT_OR_DBL  *get_gquad_pf_matrix( short *S,
                                  FLT_OR_DBL *scale,
                                  pf_paramT *pf);

void        get_gquad_pattern_mfe(short *S,
                                  int i,
                                  int j,
                                  paramT *P,
                                  int *L,
                                  int l[3]);

void        get_gquad_pattern_pf( short *S,
                                  int i,
                                  int j,
                                  pf_paramT *pf,
                                  int *L,
                                  int l[3]);

plist       *get_plist_gquad_from_pr( short *S,
                                      int gi,
                                      int gj,
                                      FLT_OR_DBL *G,
                                      FLT_OR_DBL *probs,
                                      FLT_OR_DBL *scale,
                                      pf_paramT *pf);
plist       *get_plist_gquad_from_pr_max(short *S,
                                      int gi,
                                      int gj,
                                      FLT_OR_DBL *G,
                                      FLT_OR_DBL *probs,
                                      FLT_OR_DBL *scale,
                                      int *L,
                                      int l[3],
                                      pf_paramT *pf);

plist       *get_plist_gquad_from_db( const char *structure,
                                      float pr);

/**
 *  given a dot-bracket structure (possibly) containing gquads encoded
 *  by '+' signs, find first gquad, return end position or 0 if none found
 *  Upon return L and l[] contain the number of stacked layers, as well as
 *  the lengths of the linker regions.  
 *  To parse a string with many gquads, call parse_gquad repeatedly e.g.
 *  end1 = parse_gquad(struc, &L, l); ... ;
 *  end2 = parse_gquad(struc+end1, &L, l); end2+=end1; ... ;
 *  end3 = parse_gquad(struc+end2, &L, l); end3+=end2; ... ; 
 */
int         parse_gquad(const char *struc, int *L, int l[3]);



#endif
