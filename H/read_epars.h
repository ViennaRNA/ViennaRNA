#ifndef __VIENNA_RNA_PACKAGE_READ_EPARS_H__
#define __VIENNA_RNA_PACKAGE_READ_EPARS_H__

enum parset {
  UNKNOWN= -1, QUIT,
  S, S_H, HP, HP_H, B, B_H, IL, IL_H, MMH, MMH_H, MMI, MMI_H,
  MMI1N, MMI1N_H, MMI23, MMI23_H, MMM, MMM_H, MME, MME_H, D5, D5_H, D3, D3_H,
  INT11, INT11_H, INT21, INT21_H, INT22, INT22_H, ML, TL,
  TRI, HEX, TE, NIN, MISC,DUMP, HELP}; 

void  read_parameter_file(const char fname[]);
void  write_parameter_file(const char fname[]);


#endif
