/*
   prototypes for energy_par.c
*/

#ifndef VIENNA_RNA_PACKAGE_PARAMS_DEFAULT_H
#define VIENNA_RNA_PACKAGE_PARAMS_DEFAULT_H

#include <ViennaRNA/params/constants.h>

#define PUBLIC


extern double lxc37;   /* parameter for logarithmic loop
			  energy extrapolation            */

extern int stack37[NBPAIRS+1][NBPAIRS+1];
extern int stackdH[NBPAIRS+1][NBPAIRS+1]; /* stack enthalpies */

extern int hairpin37[31];
extern int hairpindH[31];
extern int bulge37[31];
extern int bulgedH[31];
extern int internal_loop37[31];
extern int internal_loopdH[31];
extern int mismatchI37[NBPAIRS+1][5][5];  /* internal loop mismatches */
extern int mismatchIdH[NBPAIRS+1][5][5];  /* internal loop mismatches */
extern int mismatch1nI37[NBPAIRS+1][5][5];  /* internal loop mismatches */
extern int mismatch23I37[NBPAIRS+1][5][5];  /* internal loop mismatches */
extern int mismatch1nIdH[NBPAIRS+1][5][5];  /* internal loop mismatches */
extern int mismatch23IdH[NBPAIRS+1][5][5];  /* internal loop mismatches */
extern int mismatchH37[NBPAIRS+1][5][5];  /* same for hairpins */
extern int mismatchM37[NBPAIRS+1][5][5];  /* same for multiloops */
extern int mismatchHdH[NBPAIRS+1][5][5];  /* same for hairpins */
extern int mismatchMdH[NBPAIRS+1][5][5];  /* same for multiloops */
extern int mismatchExt37[NBPAIRS+1][5][5];
extern int mismatchExtdH[NBPAIRS+1][5][5];

extern int dangle5_37[NBPAIRS+1][5];      /* 5' dangle exterior of pair */
extern int dangle3_37[NBPAIRS+1][5];      /* 3' dangle */
extern int dangle3_dH[NBPAIRS+1][5];       /* corresponding enthalpies */
extern int dangle5_dH[NBPAIRS+1][5];

extern int int11_37[NBPAIRS+1][NBPAIRS+1][5][5]; /* 1x1 internal loops */
extern int int11_dH[NBPAIRS+1][NBPAIRS+1][5][5];

extern int int21_37[NBPAIRS+1][NBPAIRS+1][5][5][5]; /* 2x1 internal loops */
extern int int21_dH[NBPAIRS+1][NBPAIRS+1][5][5][5];

extern int int22_37[NBPAIRS+1][NBPAIRS+1][5][5][5][5]; /* 2x2 internal loops */
extern int int22_dH[NBPAIRS+1][NBPAIRS+1][5][5][5][5];

/* constants for linearly destabilizing contributions for multi-loops
   F = ML_closing + ML_intern*(k-1) + ML_BASE*u  */
extern int ML_BASE37;
extern int ML_BASEdH;
extern int ML_closing37;
extern int ML_closingdH;
extern int ML_intern37;
extern int ML_interndH;

extern int TripleC37;
extern int TripleCdH;
extern int MultipleCA37;
extern int MultipleCAdH;
extern int MultipleCB37;
extern int MultipleCBdH;

/* Ninio-correction for asymmetric internal loops with branches n1 and n2 */
/*    ninio_energy = min{max_ninio, |n1-n2|*F_ninio[min{4.0, n1, n2}] } */
extern int  MAX_NINIO;                   /* maximum correction */
extern int ninio37;
extern int niniodH;
/* penalty for helices terminated by AU (actually not GC) */
extern int TerminalAU37;
extern int TerminalAUdH;
/* penalty for forming bi-molecular duplex */
extern int DuplexInit37;
extern int DuplexInitdH;
/* stabilizing contribution due to special hairpins of size 4 (tetraloops) */
extern char Tetraloops[281];  /* string containing the special tetraloops */
extern int  Tetraloop37[40];  /* Bonus energy for special tetraloops */
extern int  TetraloopdH[40];
extern char Triloops[241];    /* string containing the special triloops */
extern int  Triloop37[40]; /* Bonus energy for special Triloops */
extern int  TriloopdH[40]; /* Bonus energy for special Triloops */
extern char Hexaloops[361];    /* string containing the special triloops */
extern int  Hexaloop37[40]; /* Bonus energy for special Triloops */
extern int  HexaloopdH[40]; /* Bonus energy for special Triloops */

extern int GQuadAlpha37;
extern int GQuadAlphadH;
extern int GQuadBeta37;
extern int GQuadBetadH;
extern int GQuadLayerMismatch37;  /* penalty per incompatible gquad layer in a sub-alignment (applied twice for inner layers) */
extern int GQuadLayerMismatchH;
extern int GQuadLayerMismatchMax; /* maximum number of mismatching sequences in the alignment when gquad should be formed */

extern double Tmeasure;       /* temperature of param measurements */

#endif
