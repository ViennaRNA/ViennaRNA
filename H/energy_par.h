/* 
    Current energy set taken from:
  
    D.H. Turner, N. Sugimoto, and S.M. Freier.
    "RNA Structure Prediction",
    Ann. Rev. Biophys. Biophys. Chem. 17, 167-192, 1988.

    John A.Jaeger, Douglas H.Turner, and Michael Zuker.
    "Improved predictions of secondary structures for RNA",
    PNAS, 86, 7706-7710, October 1989.

    GU stacking energies from
    L. He, R. Kierzek, J. SantaLucia, A.E. Walter, D.H. Turner
    "Nearest-Neughbor Parameters for GU Mismatches...."
    Biochemistry 1991, 30 11124-11132
    
*/

#define GASCONST 1.98717  /* in [cal/K] */
#define K0  273.15
#define INF 1000000
#define FORBIDDEN 9999
#define BONUS 10000
#define NBPAIRS 7         /* possible pairs are AU UA GC CG GU UG */
#define TURN 3            /* minimal Hairpin size */
#define MAXLOOP 30        /* maximal size for interior loops */
static double lxc37=107.856;   /* parameter for logarithmic loop
				  energy extrapolation            */

#define NST 0   /* Energy for nonstandard stacked pairs */
static int enthalpies[NBPAIRS+1][NBPAIRS+1] = 
   /* enthalpies (0.01*kcal/mol at 37 C) for stacked pairs */
   /*              CG      GC      GU      UG      AU      UA    @  */
   {  {    INF,    INF,    INF,    INF,    INF,    INF,    INF, INF}, 
      {    INF,  -1220,   -800,   -310,  -1120,  -1050,   -760, NST},
      {    INF,  -1420,  -1220,   -630,   -950,  -1330,  -1020, NST},
      {    INF,   -950,  -1120,  -2000,  -1605,  -1360,   -400, NST},
      {    INF,   -630,   -310,  -1110,  -2000,   -850,   -270, NST},
      {    INF,  -1020,   -760,   -270,   -400,   -660,   -570, NST},
      {    INF,  -1330,  -1050,   -850,  -1360,   -810,   -660, NST},
      {	   INF,    NST,    NST,    NST,    NST,    NST,    NST, NST}};

static int entropies[NBPAIRS+1][NBPAIRS+1] = 
   /* entropies (0.1*cal/(K mol)) for stacked pairs  */
   /*              CG      GC      GU      UG      AU      UA    @   */ 
   {  {   -INF,   -INF,   -INF,   -INF,   -INF,   -INF,   -INF, -INF },
      {   -INF,   -297,   -194,    -62,   -301,   -278,   -192,   0  },
      {   -INF,   -349,   -297,   -158,   -239,   -355,   -262,   0  },
      {   -INF,   -239,   -301,   -631,   -507,   -402,    -97,   0  },
      {   -INF,   -158,    -62,   -349,   -631,   -249,   - 71,   0  },
      {   -INF,   -262,   -192,   - 71,    -97,   -184,   -155,   0  },
      {   -INF,   -355,   -278,   -249,   -402,   -226,   -184,   0  },
      {   -INF,      0,      0,      0,      0,      0,      0,   0}};

static int hairpin37[31] = {
   INF, INF, INF, 450, 550, 490, 510, 520, 550, 580, 591,
        602, 611, 620, 628, 635, 642, 649, 655, 661, 666,
        671, 676, 681, 686, 690, 694, 698, 702, 706, 710};

static int bulge37[31] = {
   INF, 390, 310, 350, 420, 480, 500, 516, 531, 543, 555,
        565, 574, 583, 591, 598, 605, 612, 618, 624, 630,
        635, 640, 645, 649, 654, 658, 662, 666, 670, 673};

static int internal_loop37[31] = {
   INF, INF, 410, 450, 490, 530, 570, 587, 601, 614, 625,
        635, 645, 653, 661, 669, 676, 682, 688, 694, 700,
        705, 710, 715, 720, 724, 728, 732, 736, 740, 744};

#define DEF -50    /* Default terminal mismatch, used for I and any non_pairing bases */

/* terminal mismatches */
#define NSM 0  /* terminal mismatch for non standard pairs */
static int mismatch[5][5][5][5]=                   /*  2   */
/*             @@  @A  @C  @G  @U    A@   AA   AC   AG    AU    C@   CA   CC   CG   CU     G@   GA   GC   GG    GU    U@   UA    UC    UG    UU   */
/*  @@  */{{{{DEF,DEF,DEF,DEF,DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}},
/*  @A  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}},
/*  @C  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}},
/*  @G  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}},
/*  @U  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}}},
/*  A@  */ {{{DEF,DEF,DEF,DEF,DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}},
/*  AA  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, NSM, NSM, NSM,-100}, {DEF, NSM, NSM,-110, NSM}, {DEF, NSM,-190, NSM,-150}, {DEF, -80, NSM, -80, NSM}},
/*  AC  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, NSM, NSM, NSM, -70}, {DEF, NSM, NSM,-110, NSM}, {DEF, NSM,-100, NSM, -90}, {DEF, -70, NSM, -70, NSM}},
/*  AG  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, NSM, NSM, NSM,-110}, {DEF, NSM, NSM,-160, NSM}, {DEF, NSM,-190, NSM,-150}, {DEF, -80, NSM, -80, NSM}},
/*  AU  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, -80,-100,-100,-100}, {DEF, -70, -70, -70, -70}, {DEF, -80,-100,-100,-100}, {DEF, -80, -80, -80, -80}}},
/*  C@  */ {{{DEF,DEF,DEF,DEF,DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}},
/*  CA  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, NSM, NSM, NSM, -80}, {DEF, NSM, NSM,-130, NSM}, {DEF, NSM,-200, NSM,-140}, {DEF,-100, NSM,-100, NSM}},
/*  CC  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, NSM, NSM, NSM, -60}, {DEF, NSM, NSM, -60, NSM}, {DEF, NSM,-110, NSM, -90}, {DEF, -70, NSM, -70, NSM}},
/*  CG  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF,-190,-200,-190,-190}, {DEF,-100,-110,-100, -80}, {DEF,-190,-190,-190,-190}, {DEF,-140,-150,-140,-120}},
/*  CU  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, NSM, NSM, NSM, -60}, {DEF, NSM, NSM, -80, NSM}, {DEF, NSM,-150, NSM,-110}, {DEF, -80, NSM, -80, NSM}}},
/*  C@  */ {{{DEF,DEF,DEF,DEF,DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}},
/*  GA  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, NSM, NSM, NSM,-110}, {DEF, NSM, NSM,-130, NSM}, {DEF, NSM,-190, NSM,-150}, {DEF,-100, NSM,-100, NSM}},
/*  GC  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF,-110,-130,-130,-130}, {DEF,-110, -60, -60, -50}, {DEF,-160,-150,-140,-150}, {DEF, -80, -80, -80, -70}},
/*  GG  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, NSM, NSM, NSM,-120}, {DEF, NSM, NSM,-140, NSM}, {DEF, NSM,-190, NSM,-160}, {DEF,-100, NSM,-100, NSM}},
/*  GU  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, -80,-100,-100,-100}, {DEF,-70,  -70, -70, -70}, {DEF, -80,-100,-100, -100}, {DEF, -80, -80, -80, -80}}},
/*  U@  */ {{{DEF,DEF,DEF,DEF,DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}, {DEF, DEF, DEF, DEF, DEF}},
/*  UA  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF,-100, -80,-110, -90}, {DEF,-70,  -60, -30, -50}, {DEF,-110, -90,-120, -90}, {DEF, -30, -60, -30, -50}},
/*  UC  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, NSM, NSM, NSM, -50}, {DEF, NSM, NSM, -50, NSM}, {DEF, NSM, -80, NSM, -70}, {DEF, -70, NSM, -70, NSM}},
/*  UG  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF,-150,-140,-150,-140}, {DEF,-90,  -90, -70, -70}, {DEF,-150,-140,-160,-140}, {DEF, -90,-110, -90, -90}},
/*  UU  */  {{DEF,DEF,DEF,DEF,DEF}, {DEF, NSM, NSM, NSM, -50}, {DEF, NSM, NSM, -70, NSM}, {DEF, NSM,-120, NSM, -90}, {DEF, -80, NSM, -80, NSM}}}};

/* 5' dangling ends */
static dangle5[NBPAIRS+1][5]=
{/*    @    A    C    G    U   */
   { INF, INF, INF, INF, INF}, /* no pair */
   {   0, -50, -20, -20, -10}, /* CG   stacks on C */
   {   0, -20, -30, -00, -00}, /* GC   stacks on G */
   {   0, -20, -20, -20, -20}, /* GU */
   {   0, -20, -20, -20, -20}, /* UG */
   {   0, -30, -30, -40, -20}, /* AU */
   {   0, -30, -10, -20, -20},  /* UA */
   {   0,   0,   0,   0,   0}  /*  @ */
};

/* 3' dangling ends */
static dangle3[NBPAIRS+1][5]=
{/*    @     A     C     G     U   */
   { INF,  INF,  INF,  INF,  INF},  /* no pair */
   { DEF, -110,  -40, -130,  -60},  /* CG  stacks on G */
   { DEF, -170,  -80, -170, -120},  /* GC */      
   { DEF, -120,  -50, -120,  -70},  /* GU */
   { DEF,  -80,  -50,  -80,  -60},  /* UG */
   { DEF,  -70,  -10,  -70,  -10},  /* AU */
   { DEF,  -80,  -50,  -80,  -60},  /* UA */
   {   0,    0,     0,    0,   0}   /*  @ */
};


/* constants for linearly destabilizing contributions for multi-loops
   F = ML_closing + ML_intern*(k-1) + ML_BASE*u  */
static int ML_BASE37 = 40;
static int ML_closing37[NBPAIRS+1] = { INF, 460, 460, 460, 460, 460, 460, 460};
static int ML_intern37[NBPAIRS+1] =  { INF,  10,  10,  10,  10,  10,  10,  10};
/*                                   NoPair  CG   GC   GU   UG   AU   UA   @ */ 

/* Ninio-correction for asymmetric internal loops with branches n1 and n2 */
/*    ninio_energy = min{max_ninio, |n1-n2|*F_ninio[min{4.0, n1, n2}] } */
static int         MAX_NINIO = 300;                   /* maximum correction */
static int F_ninio37[5] = { 0, 40, 30, 20, 10 };

/* stabilizing controbution due to special hairpins of size 4 (tetraloops) */
#define      N_TETRALOOPS    8
static int   TETRA_ENERGY37 = -200;
static char  Tetraloops[N_TETRALOOPS][5] = {
                        { 'G', 'A', 'A', 'A', '\0' },
			{ 'G', 'C', 'A', 'A', '\0' },
			{ 'G', 'A', 'G', 'A', '\0' },
			{ 'G', 'U', 'G', 'A', '\0' },
			{ 'G', 'G', 'A', 'A', '\0' },
			{ 'U', 'U', 'C', 'G', '\0' },
			{ 'U', 'A', 'C', 'G', '\0' },
			{ 'G', 'C', 'G', 'A', '\0' },};










