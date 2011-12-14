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

int **annotate_gquadruplexes(short  *S);



#endif
