#ifndef __VIENNA_RNA_PACKAGE_READ_EPARS_H__
#define __VIENNA_RNA_PACKAGE_READ_EPARS_H__

/**
 *  \file read_epars.h
 *  \brief Functions to read and write energy parameter sets from/to files
 */

/**
 *  \brief
 */
enum parset {
  UNKNOWN= -1, QUIT,
  S, S_H, HP, HP_H, B, B_H, IL, IL_H, MMH, MMH_H, MMI, MMI_H,
  MMI1N, MMI1N_H, MMI23, MMI23_H, MMM, MMM_H, MME, MME_H, D5, D5_H, D3, D3_H,
  INT11, INT11_H, INT21, INT21_H, INT22, INT22_H, ML, TL,
  TRI, HEX, NIN, MISC};

/**
 *  \brief Read energy parameters from a file
 * 
 *  \param fname  The path to the file containing the energy parameters
 */
void  read_parameter_file(const char fname[]);

/**
 *  \brief Write energy parameters to a file
 * 
 *  \param fname  A filename (path) for the file where the current energy parameters will be written to
 */
void  write_parameter_file(const char fname[]);

/**
 *  \brief
 */
enum  parset gettype(const char *ident);

/**
 *  \brief
 */
char  *settype(enum parset s);

#endif
