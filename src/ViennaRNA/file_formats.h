#ifndef __VIENNA_RNA_PACKAGE_FILE_FORMATS_H__
#define __VIENNA_RNA_PACKAGE_FILE_FORMATS_H__

/**
 *  \file file_formats.h
 *  \brief Various functions dealing with file formats for RNA sequences, structures, and alignments
 */

#include <stdio.h>

#include <ViennaRNA/data_structures.h>

/**
 *  \brief Print a secondary structure as helix list
 *
 *  \param  db    The structure in dot-bracket format
 *  \param  file  The file handle used to print to (print defaults to 'stdout' if(file == NULL) )
 */
void vrna_structure_print_helix_list(const char *db, FILE *file);

/**
 *  \brief Print a secondary structure as connect table
 *
 *  \param  seq         The RNA sequence
 *  \param  db          The structure in dot-bracket format
 *  \param  energy      The free energy of the structure
 *  \param  identifier  An optional identifier for the sequence
 *  \param  file  The file handle used to print to (print defaults to 'stdout' if(file == NULL) )
 */
void vrna_structure_print_ct( const char *seq,
                              const char *db,
                              float energy,
                              const char *identifier,
                              FILE *file);

/**
 *  \brief Print a secondary structure in bpseq format
 *
 *  \param  seq         The RNA sequence
 *  \param  db          The structure in dot-bracket format
 *  \param  file  The file handle used to print to (print defaults to 'stdout' if(file == NULL) )
 */
void vrna_structure_print_bpseq(const char *seq,
                                const char *db,
                                FILE *file);

#endif
