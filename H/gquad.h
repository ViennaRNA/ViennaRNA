#ifndef __VIENNA_RNA_PACKAGE_GQUAD_H__
#define __VIENNA_RNA_PACKAGE_GQUAD_H__

/**
 *  \file quad.h
 *  \brief Various functions related to G-quadruplex computations
 */

#define   VRNA_GQUAD_MAX_STACK_SIZE     7
#define   VRNA_GQUAD_MIN_STACK_SIZE     2
#define   VRNA_GQUAD_MAX_LINKER_LENGTH  25
#define   VRNA_GQUAD_MIN_LINKER_LENGTH  1
#define   VRNA_GQUAD_MIN_BOX_SIZE       (4*VRNA_GQUAD_MIN_STACK_SIZE)+(3*VRNA_GQUAD_MIN_LINKER_LENGTH)
#define   VRNA_GQUAD_MAX_BOX_SIZE       (4*VRNA_GQUAD_MAX_STACK_SIZE)+(3*VRNA_GQUAD_MAX_LINKER_LENGTH)

int **annotate_gquadruplexes(short  *S);

int *get_gquad_matrix(short *S);

int gquad_contribution(int L, int l1, int l2, int l3);

#endif
