/* 
    Current energy set taken from:
  
    D.H. Turner, N. Sugimoto, and S.M. Freier.
    "RNA Structure Prediction",
    Ann. Rev. Biophys. Biophys. Chem. 17, 167-192, 1988.

    John A.Jaeger, Douglas H.Turner, and Michael Zuker.
    "Improved predictions of secondary structures for RNA",
    PNAS, 86, 7706-7710, October 1989.
    
    GU stacking energies from:
    L. He, R. Kierzek, J. SantaLucia, A.E. Walter, D.H. Turner
    "Nearest-Neughbor Parameters for GU Mismatches...."
    Biochemistry 1991, 30 11124-11132

    asymetry corrections:
    A.E. Peritz, R. Kierzek, N, Sugimoto, D.H. Turner
    "Thermodynamic Study of Internal Loops in Oligoribonucleotides..."
    Biochemistry 1991, 30, 6428--6435

    See also:
    A. Walter, D Turner, J Kim, M Lyttle, P M"uller, D Mathews, M Zuker
    "Coaxial stckaing of helices enhances binding of oligoribonucleotides.."
    PNAS, 91, pp 9218-9222, 1994

    Unpublished parameters taken from mfold-2.3
*/

#include "energy_const.h"
static char rcsid[] = "$Id: energy_par.c,v 1.1 1997/10/27 13:22:33 ivo Exp $";

#define NST 0     /* Energy for nonstandard stacked pairs */
#define DEF -50   /* Default terminal mismatch, used for I */
                  /* and any non_pairing bases */
#define NSM 0     /* terminal mismatch for non standard pairs */
 
#define PUBLIC

PUBLIC double Tmeasure = 37+K0;  /* temperature of param measurements */
PUBLIC double lxc37=107.856;     /* parameter for logarithmic loop
				    energy extrapolation            */

PUBLIC int stack37[NBPAIRS+1][NBPAIRS+1] =
/*          CG     GC     GU     UG     AU     UA  */
{ {  INF,   INF,   INF,   INF,   INF,   INF,   INF, INF},
  {  INF,  -290,  -200,  -120,  -190,  -180,  -170, NST},
  {  INF,  -340,  -290,  -140,  -210,  -230,  -210, NST},
  {  INF,  -210,  -190,   -40,   150,  -110,  -100, NST},
  {  INF,  -140,  -120,   -20,   -40,   -80,   -50, NST},
  {  INF,  -210,  -170,   -50,  -100,   -90,   -90, NST},
  {  INF,  -230,  -180,   -80,  -110,  -110,   -90, NST},
  {  INF,   NST,   NST,   NST,   NST,   NST,   NST, NST}};

/* enthalpies (0.01*kcal/mol at 37 C) for stacked pairs */
/* different from mfold-2.3, which values from mfold-2.2 */
PUBLIC int enthalpies[NBPAIRS+1][NBPAIRS+1] = 
/*          CG     GC     GU     UG     AU     UA  */
{ {  INF,   INF,   INF,   INF,   INF,   INF,   INF, INF}, 
  {  INF, -1220,  -800,  -310, -1120, -1050,  -760, NST},
  {  INF, -1420, -1220,  -630,  -950, -1330, -1020, NST},
  {  INF,  -950, -1120, -2000, -1605, -1360,  -400, NST},
  {  INF,  -630,  -310, -1110, -2000,  -850,  -270, NST},
  {  INF, -1020,  -760,  -270,  -400,  -660,  -570, NST},
  {  INF, -1330, -1050,  -850, -1360,  -810,  -660, NST},
  {  INF,   NST,   NST,   NST,   NST,   NST,   NST, NST}};


PUBLIC int entropies[NBPAIRS+1][NBPAIRS+1] = /* not used anymore */
   /* entropies (0.1*cal/(K mol)) for stacked pairs  */
   /*           CG     GC     GU     UG     AU     UA  */ 
   {  {   -INF,   -INF,   -INF,   -INF,   -INF,   -INF,   -INF, -INF },
      {   -INF,   -297,   -194,    -62,   -301,   -278,   -192,   0  },
      {   -INF,   -349,   -297,   -158,   -239,   -355,   -262,   0  },
      {   -INF,   -239,   -301,   -631,   -507,   -402,    -97,   0  },
      {   -INF,   -158,    -62,   -349,   -631,   -249,   - 71,   0  },
      {   -INF,   -262,   -192,   - 71,    -97,   -184,   -155,   0  },
      {   -INF,   -355,   -278,   -249,   -402,   -226,   -184,   0  },
      {   -INF,      0,      0,      0,      0,      0,      0,   0}};

/* old values are here just for comparison */
PUBLIC int oldhairpin37[31] = {
   INF, INF, INF, 450, 550, 490, 510, 520, 550, 580, 591,
        602, 611, 620, 628, 635, 642, 649, 655, 661, 666,
        671, 676, 681, 686, 690, 694, 698, 702, 706, 710};

PUBLIC int hairpin37[31] = {
   INF, INF, INF, 410, 490, 440, 470, 500, 510, 520, 531,
        542, 551, 560, 568, 575, 582, 589, 595, 601, 606,
        611, 616, 621, 626, 630, 634, 638, 642, 646, 650};

PUBLIC int bulge37[31] = {
   INF, 390, 310, 350, 420, 480, 500, 516, 531, 543, 555,
        565, 574, 583, 591, 598, 605, 612, 618, 624, 630,
        635, 640, 645, 649, 654, 658, 662, 666, 670, 673};

/* internal_loop[3] changed from 4.5 to 5.1 */
PUBLIC int internal_loop37[31] = {
   INF, INF, 410, 510, 490, 530, 570, 587, 601, 614, 625,
        635, 645, 653, 661, 669, 676, 682, 688, 694, 700,
        705, 710, 715, 720, 724, 728, 732, 736, 740, 744};


/* terminal mismatches */
/* internal loops, free energies at 37C */
PUBLIC int old_mismatch_37[NBPAIRS+1][5][5] =
{ /* CG */
  {{  0,   0,    0,   0,   0}, /* @@  @A  @C  @G  @U */
   {DEF, -190,-200,-190,-190}, /* AA   AC   AG   AU  A@ */
   {DEF, -100,-110,-100, -80}, /* CA   CC   CG   CU  C@ */
   {DEF, -190,-190,-190,-190}, /* GA   GC   GG   GU  G@ */
   {DEF, -140,-150,-140,-120}, /* UA   UC   UG   UU  U@ */
  },
  /* GC */
  {{  0,   0,    0,   0,   0}, /* @@  @A  @C  @G  @U */
   {DEF, -110,-130,-130,-130}, /* AA   AC   AG   AU  A@ */
   {DEF, -110, -60, -60, -50}, /* CA   CC   CG   CU  C@ */
   {DEF, -160,-150,-140,-150}, /* GA   GC   GG   GU  G@ */
   {DEF,  -80, -80, -80, -70}, /* UA   UC   UG   UU  U@ */
  }, 
  /* GU */
  {{  0,   0,    0,   0,   0}, /* @@  @A  @C  @G  @U */
   {DEF,  -80,-100,-100,-100}, /* AA   AC   AG   AU  A@ */
   {DEF,  -70, -70, -70, -70}, /* CA   CC   CG   CU  C@ */
   {DEF,  -80,-100,-100,-100}, /* GA   GC   GG   GU  G@ */
   {DEF,  -80, -80, -80, -80}, /* UA   UC   UG   UU  U@ */
  },
  /* UG */
  {{  0,   0,    0,   0,   0}, /* @@  @A  @C  @G  @U */
   {DEF, -150,-140,-150,-140}, /* AA   AC   AG   AU  A@ */
   {DEF,  -90, -90, -70, -70}, /* CA   CC   CG   CU  C@ */
   {DEF, -150,-140,-160,-140}, /* GA   GC   GG   GU  G@ */
   {DEF,  -90,-110, -90, -90}, /* UA   UC   UG   UU  U@ */
  },
  /* AU */
  {{  0,   0,    0,   0,   0}, /* @@  @A  @C  @G  @U */
   {DEF,  -80,-100,-100,-100}, /* AA   AC   AG   AU  A@ */
   {DEF,  -70, -70, -70, -70}, /* CA   CC   CG   CU  C@ */
   {DEF,  -80,-100,-100,-100}, /* GA   GC   GG   GU  G@ */
   {DEF,  -80, -80, -80, -80}, /* UA   UC   UG   UU  U@ */
  },
  /* UA */
  {{  0,   0,    0,   0,   0}, /* @@  @A  @C  @G  @U */
   {DEF, -100, -80,-110, -90}, /* AA   AC   AG   AU  A@ */
   {DEF,  -70, -60, -30, -50}, /* CA   CC   CG   CU  C@ */
   {DEF, -110, -90,-120, -90}, /* GA   GC   GG   GU  G@ */
   {DEF,  -30, -60, -30, -50}, /* UA   UC   UG   UU  U@ */
  }
};

/* mismatch free energies for interior loops at 37C */
PUBLIC int mismatchI37[NBPAIRS+1][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -150, -150, -270, -190}, /* A@  AA  AC  AG  AU */
   { DEF, -150, -150, -100, -150}, /* C@  CA  CC  CG  CU */
   { DEF, -270, -190, -150, -190}, /* G@  GA  GC  GG  GU */
   { DEF, -140, -150, -140, -250}},/* U@  UA  UC  UG  UU */
  { /* GC */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -150, -150, -270, -130}, /* A@  AA  AC  AG  AU */
   { DEF, -150, -150,  -60, -150}, /* C@  CA  CC  CG  CU */
   { DEF, -270, -150, -150, -150}, /* G@  GA  GC  GG  GU */
   { DEF,  -80, -150,  -80, -250}},/* U@  UA  UC  UG  UU */
  { /* GU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -150, -150, -270, -130}, /* A@  AA  AC  AG  AU */
   { DEF, -150, -150,  -60, -150}, /* C@  CA  CC  CG  CU */
   { DEF, -270, -150, -150, -150}, /* G@  GA  GC  GG  GU */
   { DEF,  -80, -150,  -80, -250}},/* U@  UA  UC  UG  UU */
  { /* UG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -100, -100, -220,  -40}, /* A@  AA  AC  AG  AU */
   { DEF, -100, -100,   20, -100}, /* C@  CA  CC  CG  CU */
   { DEF, -220,  -40, -100,  -40}, /* G@  GA  GC  GG  GU */
   { DEF,   20, -100,   20, -200}},/* U@  UA  UC  UG  UU */
  { /* AU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -100, -100, -220,  -50}, /* A@  AA  AC  AG  AU */
   { DEF, -100, -100,  -20, -100}, /* C@  CA  CC  CG  CU */
   { DEF, -220,  -50, -100,  -50}, /* G@  GA  GC  GG  GU */
   { DEF,  -30, -100,  -30, -200}},/* U@  UA  UC  UG  UU */
  { /* UA */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -100, -100, -220,  -40}, /* A@  AA  AC  AG  AU */
   { DEF, -100, -100,   20, -100}, /* C@  CA  CC  CG  CU */
   { DEF, -220,  -40, -100,  -40}, /* G@  GA  GC  GG  GU */
   { DEF,   20, -100,   20, -200}},/* U@  UA  UC  UG  UU */
  { /* @@ */
   {DEF,DEF,DEF,DEF,DEF},{DEF,DEF,DEF,DEF,DEF},{DEF,DEF,DEF,DEF,DEF},
   {DEF,DEF,DEF,DEF,DEF},{DEF,DEF,DEF,DEF,DEF}}
};

/* mismatch free energies for hairpins at 37C */
PUBLIC int mismatchH37[NBPAIRS+1][5][5] =
{ /* @@ */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -140, -200, -210, -190}, /* A@  AA  AC  AG  AU */
   { DEF, -100, -110, -100,  -80}, /* C@  CA  CC  CG  CU */
   { DEF, -210, -190, -140, -190}, /* G@  GA  GC  GG  GU */
   { DEF, -140, -150, -140, -120}},/* U@  UA  UC  UG  UU */
  { /* GC */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -110, -130, -200, -130}, /* A@  AA  AC  AG  AU */
   { DEF, -110,  -60,  -60,  -50}, /* C@  CA  CC  CG  CU */
   { DEF, -230, -150, -140, -150}, /* G@  GA  GC  GG  GU */
   { DEF,  -80,  -80,  -80,  -70}},/* U@  UA  UC  UG  UU */
  { /* GU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF,  -80, -100, -170, -100}, /* A@  AA  AC  AG  AU */
   { DEF,  -70,  -70,  -70,  -70}, /* C@  CA  CC  CG  CU */
   { DEF, -150, -100, -100, -100}, /* G@  GA  GC  GG  GU */
   { DEF,  -80,  -80,  -80,  -80}},/* U@  UA  UC  UG  UU */
  { /* UG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -120, -140, -200, -140}, /* A@  AA  AC  AG  AU */
   { DEF,  -90,  -90,  -70,  -70}, /* C@  CA  CC  CG  CU */
   { DEF, -200, -140, -130, -140}, /* G@  GA  GC  GG  GU */
   { DEF,  -90, -110,  -90,  -90}},/* U@  UA  UC  UG  UU */
  { /* AU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF,  -80, -100, -170, -100}, /* A@  AA  AC  AG  AU */
   { DEF,  -70,  -70,  -70,  -70}, /* C@  CA  CC  CG  CU */
   { DEF, -150, -100, -100, -100}, /* G@  GA  GC  GG  GU */
   { DEF,  -80,  -80,  -80,  -80}},/* U@  UA  UC  UG  UU */
  { /* UA */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -100,  -80, -180,  -90}, /* A@  AA  AC  AG  AU */
   { DEF,  -70,  -60,  -30,  -50}, /* C@  CA  CC  CG  CU */
   { DEF, -180,  -90, -120,  -90}, /* G@  GA  GC  GG  GU */
   { DEF,  -30,  -60,  -30,  -50}},/* U@  UA  UC  UG  UU */
  { /* @@ */
   {DEF,DEF,DEF,DEF,DEF},{DEF,DEF,DEF,DEF,DEF},{DEF,DEF,DEF,DEF,DEF},
   {DEF,DEF,DEF,DEF,DEF},{DEF,DEF,DEF,DEF,DEF}}
};

/* mismatch enthalpies for temperature scaling */
PUBLIC int mism_H[NBPAIRS+1][5][5] =
{ /* no pair */
  {{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0},{0,0,0,0,0}},
  { /* CG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF,-1030, -950,-1030,-1030}, /* A@  AA  AC  AG  AU */
   { DEF, -520, -450, -520, -670}, /* C@  CA  CC  CG  CU */
   { DEF, -940, -940, -940, -940}, /* G@  GA  GC  GG  GU */
   { DEF, -810, -740, -810, -860}},/* U@  UA  UC  UG  UU */
  { /* GC */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -520, -880, -560, -880}, /* A@  AA  AC  AG  AU */
   { DEF, -720, -310, -310, -390}, /* C@  CA  CC  CG  CU */
   { DEF, -710, -740, -620, -740}, /* G@  GA  GC  GG  GU */
   { DEF, -500, -500, -500, -570}},/* U@  UA  UC  UG  UU */
  { /* GU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -430, -600, -600, -600}, /* A@  AA  AC  AG  AU */
   { DEF, -260, -240, -240, -240}, /* C@  CA  CC  CG  CU */
   { DEF, -340, -690, -690, -690}, /* G@  GA  GC  GG  GU */
   { DEF, -330, -330, -330, -330}},/* U@  UA  UC  UG  UU */
  { /* UG */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -720, -790, -960, -810}, /* A@  AA  AC  AG  AU */
   { DEF, -480, -480, -360, -480}, /* C@  CA  CC  CG  CU */
   { DEF, -660, -810, -920, -810}, /* G@  GA  GC  GG  GU */
   { DEF, -550, -440, -550, -360}},/* U@  UA  UC  UG  UU */
  { /* AU */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -430, -600, -600, -600}, /* A@  AA  AC  AG  AU */
   { DEF, -260, -240, -240, -240}, /* C@  CA  CC  CG  CU */
   { DEF, -340, -690, -690, -690}, /* G@  GA  GC  GG  GU */
   { DEF, -330, -330, -330, -330}},/* U@  UA  UC  UG  UU */
  { /* UA */
   {   0,    0,    0,    0,    0}, /* @@  @A  @C  @G  @U */
   { DEF, -400, -630, -890, -590}, /* A@  AA  AC  AG  AU */
   { DEF, -430, -510, -200, -180}, /* C@  CA  CC  CG  CU */
   { DEF, -380, -680, -890, -680}, /* G@  GA  GC  GG  GU */
   { DEF, -280, -140, -280, -140}},/* U@  UA  UC  UG  UU */
  { /* nonstandard pair */
   {DEF,DEF,DEF,DEF,DEF},{DEF,DEF,DEF,DEF,DEF},{DEF,DEF,DEF,DEF,DEF},
   {DEF,DEF,DEF,DEF,DEF},{DEF,DEF,DEF,DEF,DEF}}
};

/* 5' dangling ends (dangle is exterior of base pair) */
PUBLIC int dangle5_37[NBPAIRS+1][5]=
{/*    @    A    C    G    U   */
   { INF, INF, INF, INF, INF}, /* no pair */
   {   0, -50, -30, -20, -10}, /* CG   stacks on C */
   {   0, -20, -30,   0,   0}, /* GC   stacks on G */
   {   0, -20, -20, -20, -20}, /* GU */
   {   0, -20, -20, -20, -20}, /* UG */
   {   0, -30, -30, -40, -20}, /* AU */
   {   0, -30, -10, -20, -20}, /* UA */
   {   0,   0,   0,   0,   0}  /*  @ */
};

/* 3' dangling ends */
PUBLIC int dangle3_37[NBPAIRS+1][5]=
{/*    @     A     C     G     U   */
   { INF,  INF,  INF,  INF,  INF},  /* no pair */
   { DEF, -110,  -40, -130,  -60},  /* CG  stacks on G */
   { DEF, -170,  -80, -170, -120},  /* GC */      
   { DEF, -120,  -50, -120,  -70},  /* GU */
   { DEF,  -80,  -50,  -80,  -60},  /* UG */
   { DEF,  -70,  -10,  -70,  -10},  /* AU */
   { DEF,  -80,  -50,  -80,  -60},  /* UA */
   {   0,    0,    0,    0,   0}    /*  @ */
};
				       
/* enthalpies for temperature scaling */
PUBLIC int dangle3_H[NBPAIRS+1][5] =
{/*   @     A     C     G     U   */
   { INF,  INF,  INF,  INF,  INF},  /* no pair */
   {   0, -740, -280, -640, -360},
   {   0, -900, -410, -860, -750},
   {   0, -740, -240, -720, -490},
   {   0, -490,  -90, -550, -230},
   {   0, -570,  -70, -580, -220},
   {   0, -490,  -90, -550, -230},
   {   0,    0,    0,    0,   0}
};

PUBLIC int dangle5_H[NBPAIRS+1][5] =
{/*   @     A     C     G     U   */
   { INF,  INF,  INF,  INF,  INF},  /* no pair */
   {   0, -240,  330,   80, -140},
   {   0, -160,   70, -460,  -40},
   {   0,  160,  220,   70,  310},
   {   0, -150,  510,   10,  100},
   {   0,  160,  220,   70,  310},
   {   0,  -50,  690,  -60,  -60},
   {   0,    0,    0,    0,   0}
};


/* constants for linearly destabilizing contributions for multi-loops
   F = ML_closing + ML_intern*k + ML_BASE*u  */
/* old versions erroneously used ML_intern*(k-1) */
PUBLIC int ML_BASE37 = 40;
PUBLIC int ML_closing37 = 460;
PUBLIC int ML_intern37 =  10;

/* Ninio-correction for asymmetric internal loops with branches n1 and n2 */
/*    ninio_energy = min{max_ninio, |n1-n2|*F_ninio[min{4.0, n1, n2}] } */
PUBLIC int         MAX_NINIO = 300;                   /* maximum correction */
PUBLIC int F_ninio37[5] = { 0, 40, 30, 20, 10 };      /* not used anymore */

/* stabilizing contribution due to special hairpins of size 4 (tetraloops) */

PUBLIC char Tetraloops[200] =
"GAAA "
"GCAA "
"GAGA "
"GUGA "
"GGAA "
"UUCG "
"UACG "
"GCGA "
"UCCG "
"GUAA "
"CUUG "
"AUUU "
"UUUA ";

PUBLIC int   TETRA_ENERGY37 = -200;
PUBLIC int   TETRA_ENTH37   = -400;
