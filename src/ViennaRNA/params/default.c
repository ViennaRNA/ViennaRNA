

/*
    Automatically generated using the TurnerParser
    TurnerParser (c) 2008,2009,2010
      Christian Hoener zu Siederdissen, TBI Vienna
      choener (at) tbi.univie.ac.at

    The library enabling this can be found at:
    http://hackage.haskell.org/package/BiobaseVienna
    the program can be found at:
    (sorry, not yet)
    install using cabal: cabal install (sorry, not yet)
*/

/*
     Current free energy parameters are summarized in:

     D.H.Mathews, J. Sabina, M. ZUker, D.H. Turner
     "Expanded sequence dependence of thermodynamic parameters improves
     prediction of RNA secondary structure"
     JMB, 288, pp 911-940, 1999

     Enthalpies taken from:

     A. Walter, D Turner, J Kim, M Lyttle, P M"uller, D Mathews, M Zuker
     "Coaxial stacking of helices enhances binding of oligoribonucleotides.."
     PNAS, 91, pp 9218-9222, 1994

     D.H. Turner, N. Sugimoto, and S.M. Freier.
     "RNA Structure Prediction",
     Ann. Rev. Biophys. Biophys. Chem. 17, 167-192, 1988.

     John A.Jaeger, Douglas H.Turner, and Michael Zuker.
     "Improved predictions of secondary structures for RNA",
     PNAS, 86, 7706-7710, October 1989.

     L. He, R. Kierzek, J. SantaLucia, A.E. Walter, D.H. Turner
     "Nearest-Neighbor Parameters for GU Mismatches...."
     Biochemistry 1991, 30 11124-11132

     A.E. Peritz, R. Kierzek, N, Sugimoto, D.H. Turner
     "Thermodynamic Study of Internal Loops in Oligoribonucleotides..."
     Biochemistry 1991, 30, 6428--6435
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ViennaRNA/params/default.h"

#define NST 0     /* Energy for nonstandard stacked pairs */
#define DEF -50   /* Default terminal mismatch, used for I */
                  /* and any non_pairing bases */
#define NSM 0     /* terminal mismatch for non standard pairs */

#define PUBLIC

PUBLIC double Tmeasure = 37+K0;  /* temperature of param measurements */

/* PUBLIC double lxc37=107.9; */
PUBLIC double lxc37=107.856;
PUBLIC int ML_intern37=-90;
PUBLIC int ML_interndH=-220;
PUBLIC int ML_closing37=930;
PUBLIC int ML_closingdH=3000;
PUBLIC int ML_BASE37=0;
PUBLIC int ML_BASEdH=0;
PUBLIC int MAX_NINIO=300;
PUBLIC int ninio37=60;
PUBLIC int niniodH=320;
PUBLIC int TerminalAU37=50;
PUBLIC int TerminalAUdH=370;
PUBLIC int DuplexInit37=410;
PUBLIC int DuplexInitdH=360;
PUBLIC int TripleC37=100;
PUBLIC int TripleCdH=1860;
PUBLIC int MultipleCA37=30;
PUBLIC int MultipleCAdH=340;
PUBLIC int MultipleCB37=160;
PUBLIC int MultipleCBdH=760;

PUBLIC int GQuadAlpha37 = -1800;
PUBLIC int GQuadAlphadH = -11934;
PUBLIC int GQuadBeta37 = 1200;
PUBLIC int GQuadBetadH = 0;
PUBLIC int GQuadLayerMismatch37   = 300;
PUBLIC int GQuadLayerMismatchH    = 0;
PUBLIC int GQuadLayerMismatchMax  = 1;

PUBLIC int stack37[NBPAIRS+1][NBPAIRS+1] =
{{   INF,   INF,   INF,   INF,   INF,   INF,   INF,   INF}
,{   INF,  -240,  -330,  -210,  -140,  -210,  -210,  -140}
,{   INF,  -330,  -340,  -250,  -150,  -220,  -240,  -150}
,{   INF,  -210,  -250,   130,   -50,  -140,  -130,   130}
,{   INF,  -140,  -150,   -50,    30,   -60,  -100,    30}
,{   INF,  -210,  -220,  -140,   -60,  -110,   -90,   -60}
,{   INF,  -210,  -240,  -130,  -100,   -90,  -130,   -90}
,{   INF,  -140,  -150,   130,    30,   -60,   -90,   130}};
PUBLIC int stackdH[NBPAIRS+1][NBPAIRS+1] =
{{   INF,   INF,   INF,   INF,   INF,   INF,   INF,   INF}
,{   INF, -1060, -1340, -1210,  -560, -1050, -1040,  -560}
,{   INF, -1340, -1490, -1260,  -830, -1140, -1240,  -830}
,{   INF, -1210, -1260, -1460, -1350,  -880, -1280,  -880}
,{   INF,  -560,  -830, -1350,  -930,  -320,  -700,  -320}
,{   INF, -1050, -1140,  -880,  -320,  -940,  -680,  -320}
,{   INF, -1040, -1240, -1280,  -700,  -680,  -770,  -680}
,{   INF,  -560,  -830,  -880,  -320,  -320,  -680,  -320}};

PUBLIC int hairpin37[31] = {   INF,   INF,   INF,   540,   560,   570,   540,   600,   550,   640,   650,   660,   670,   680,   690,   690,   700,   710,   710,   720,   720,   730,   730,   740,   740,   750,   750,   750,   760,   760,   770};
PUBLIC int hairpindH[31] = {   INF,   INF,   INF,   130,   480,   360,  -290,   130,  -290,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500,   500};
PUBLIC int bulge37[31] = {   INF,   380,   280,   320,   360,   400,   440,   460,   470,   480,   490,   500,   510,   520,   530,   540,   540,   550,   550,   560,   570,   570,   580,   580,   580,   590,   590,   600,   600,   600,   610};
PUBLIC int bulgedH[31] = {   INF,  1060,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710,   710};
PUBLIC int internal_loop37[31] = {   INF,   INF,   100,   100,   110,   200,   200,   210,   230,   240,   250,   260,   270,   280,   290,   290,   300,   310,   310,   320,   330,   330,   340,   340,   350,   350,   350,   360,   360,   370,   370};
PUBLIC int internal_loopdH[31] = {   INF,   INF,   -720,   -720,  -720,  -680,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130,  -130};

PUBLIC int mismatchI37[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,   -80,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -100,     0,  -100,     0}
 ,{     0,     0,     0,     0,   -60}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,   -80,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -100,     0,  -100,     0}
 ,{     0,     0,     0,     0,   -60}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,   -10,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -30,    70,   -30,    70}
 ,{    70,    70,    70,    70,    10}
 }};
PUBLIC int mismatchIdH[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{   280,     0,     0,   280,     0}
 ,{     0,     0,     0,  -340,     0}
 ,{     0,     0,     0,     0,     0}
 ,{   280,  -760,     0,   280,     0}
 ,{     0,     0,     0,     0,  -580}
 }
,{{   280,     0,     0,   280,     0}
 ,{     0,     0,     0,  -340,     0}
 ,{     0,     0,     0,     0,     0}
 ,{   280,  -760,     0,   280,     0}
 ,{     0,     0,     0,     0,  -580}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }
,{{   790,   500,   500,   790,   500}
 ,{   500,   500,   500,   170,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   790,  -260,   500,   790,   500}
 ,{   500,   500,   500,   500,   -80}
 }};

PUBLIC int mismatchH37[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{   -80,  -100,  -110,  -100,   -80}
 ,{  -140,  -150,  -150,  -140,  -150}
 ,{   -80,  -100,  -110,  -100,   -80}
 ,{  -150,  -230,  -150,  -240,  -150}
 ,{  -100,  -100,  -140,  -100,  -210}
 }
,{{   -50,  -110,   -70,  -110,   -50}
 ,{  -110,  -110,  -150,  -130,  -150}
 ,{   -50,  -110,   -70,  -110,   -50}
 ,{  -150,  -250,  -150,  -220,  -150}
 ,{  -100,  -110,  -100,  -110,  -160}
 }
,{{    20,    20,   -20,   -10,   -20}
 ,{    20,    20,   -50,   -30,   -50}
 ,{   -10,   -10,   -20,   -10,   -20}
 ,{   -50,  -100,   -50,  -110,   -50}
 ,{   -10,   -10,   -30,   -10,  -100}
 }
,{{     0,   -20,   -10,   -20,     0}
 ,{   -30,   -50,   -30,   -60,   -30}
 ,{     0,   -20,   -10,   -20,     0}
 ,{   -30,   -90,   -30,  -110,   -30}
 ,{   -10,   -20,   -10,   -20,   -90}
 }
,{{   -10,   -10,   -20,   -10,   -20}
 ,{   -30,   -30,   -50,   -30,   -50}
 ,{   -10,   -10,   -20,   -10,   -20}
 ,{   -50,  -120,   -50,  -110,   -50}
 ,{   -10,   -10,   -30,   -10,  -120}
 }
,{{     0,   -20,   -10,   -20,     0}
 ,{   -30,   -50,   -30,   -50,   -30}
 ,{     0,   -20,   -10,   -20,     0}
 ,{   -30,  -150,   -30,  -150,   -30}
 ,{   -10,   -20,   -10,   -20,   -90}
 }
,{{    20,    20,   -10,   -10,     0}
 ,{    20,    20,   -30,   -30,   -30}
 ,{     0,   -10,   -10,   -10,     0}
 ,{   -30,   -90,   -30,  -110,   -30}
 ,{   -10,   -10,   -10,   -10,   -90}
 }};
PUBLIC int mismatchHdH[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{   560,  -570,   560,  -560,  -270}
 ,{  -560,  -910,  -560,  -560,  -560}
 ,{  -270,  -570,  -340,  -570,  -270}
 ,{   560, -1400,   560,  -920,  -560}
 ,{  -530,  -570,  -530,  -570, -1440}
 }
,{{    50,  -520,    50,  -560,  -400}
 ,{  -400,  -520,  -400,  -560,  -400}
 ,{    50,  -720,    50,  -720,  -420}
 ,{  -400, -1290,  -400,  -620,  -400}
 ,{   -30,  -720,   -30,  -720, -1080}
 }
,{{   970,   140,   970,   140,   570}
 ,{   570,    30,   570,    20,   570}
 ,{   970,   140,   970,   140,   340}
 ,{   570,  -270,   570,    20,   570}
 ,{   830,   140,   830,   140,   -50}
 }
,{{   230,   100,   230,   220,   190}
 ,{  -110,  -110,  -260,  -520,  -260}
 ,{   190,   -60,  -140,   -60,   190}
 ,{   220,   100,  -260,   220,  -260}
 ,{   230,   -60,   230,   -60,   -70}
 }
,{{   970,   140,   970,   140,   570}
 ,{   570,   -20,   570,    20,   570}
 ,{   970,   140,   970,   140,   340}
 ,{   570,  -520,   570,    20,   570}
 ,{   830,   140,   830,   140,  -380}
 }
,{{   230,   -30,   230,   -60,   190}
 ,{   -30,   -30,  -260,  -520,  -260}
 ,{   190,   -60,  -140,   -60,   190}
 ,{  -260,  -590,  -260,  -520,  -260}
 ,{   230,   -60,   230,   -60,   -70}
 }
,{{   970,   140,   970,   220,   570}
 ,{   570,    30,   570,    20,   570}
 ,{   970,   140,   970,   140,   340}
 ,{   570,   100,   570,   220,   570}
 ,{   830,   140,   830,   140,   -50}
 }};

/* mismatch_multi */
PUBLIC int mismatchM37[NBPAIRS+1][5][5] =
{{ /* NP.. */
  {   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 },
 { /* CG.. */
  {   -50,  -110,   -50,  -140,   -70}
 ,{  -110,  -110,  -110,  -160,  -110}
 ,{   -70,  -150,   -70,  -150,  -100}
 ,{  -110,  -130,  -110,  -140,  -110}
 ,{   -50,  -150,   -50,  -150,   -70}
 },
 { /* GC.. */
  {   -80,  -140,   -80,  -140,  -100}
 ,{  -100,  -150,  -100,  -140,  -100}
 ,{  -110,  -150,  -110,  -150,  -140}
 ,{  -100,  -140,  -100,  -160,  -100}
 ,{   -80,  -150,   -80,  -150,  -120}
 },
 { /* GU.. */
  {   -50,   -80,   -50,   -50,   -50}
 ,{   -50,  -100,   -70,   -50,   -70}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -110,   -70,   -80,   -70}
 ,{   -50,   -80,   -50,   -80,   -50}
 },
 { /* UG.. */
  {   -30,   -30,   -60,   -60,   -60}
 ,{   -30,   -30,   -60,   -60,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -60,  -100,   -70,  -100,   -60}
 },
 { /* AU.. */
  {   -50,   -80,   -50,   -80,   -50}
 ,{   -70,  -100,   -70,  -110,   -70}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -110,   -70,  -120,   -70}
 ,{   -50,   -80,   -50,   -80,   -50}
 },
 { /* UA.. */
  {   -60,   -80,   -60,   -80,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 },
 { /* NN.. */
  {   -30,   -30,   -50,   -50,   -50}
 ,{   -30,   -30,   -60,   -50,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -50,   -80,   -50,   -80,   -50}
 }};

/* mismatch_multi_enthalpies */
PUBLIC int mismatchMdH[NBPAIRS+1][5][5] =
{{ /* NP.. */
  {   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 },
 { /* CG.. */
  {    50,  -400,    50,  -400,   -30}
 ,{  -520,  -520,  -720,  -710,  -720}
 ,{    50,  -400,    50,  -400,   -30}
 ,{  -560,  -560,  -720,  -620,  -720}
 ,{  -400,  -400,  -420,  -400,  -500}
 },
 { /* GC.. */
  {  -270,  -560,  -270,  -560,  -530}
 ,{  -570,  -910,  -570,  -820,  -570}
 ,{  -340,  -560,  -340,  -560,  -530}
 ,{  -560,  -560,  -570,  -920,  -570}
 ,{  -270,  -560,  -270,  -560,  -860}
 },
 { /* GU.. */
  {   310,  -480,  -180,   310,   140}
 ,{   310,  -480,  -430,   310,  -430}
 ,{  -140,  -630,  -510,  -630,  -140}
 ,{  -150,  -890,  -430,  -150,  -430}
 ,{   140,  -630,  -180,  -630,   140}
 },
 { /* UG.. */
  {   600,   200,   600,   200,   460}
 ,{   -60,  -340,  -230,   -60,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -230,  -350,  -230,  -350,  -230}
 ,{   200,   200,   -30,   200,   160}
 },
 { /* AU.. */
  {   140,  -400,  -180,  -380,   140}
 ,{  -380,  -400,  -430,  -380,  -430}
 ,{  -140,  -630,  -510,  -630,  -140}
 ,{  -430,  -890,  -430,  -890,  -430}
 ,{   140,  -630,  -180,  -630,   140}
 },
 { /* UA.. */
  {   600,   200,   600,   200,   460}
 ,{  -230,  -390,  -230,  -310,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -230,  -350,  -230,  -350,  -230}
 ,{   200,   200,   -30,   200,  -170}
 },
 { /* NN.. */
  {   600,   200,   600,   310,   460}
 ,{   310,  -340,  -230,   310,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -150,  -350,  -230,  -150,  -230}
 ,{   200,   200,   -30,   200,   160}
 }};

PUBLIC int mismatch1nI37[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 }};
PUBLIC int mismatch1nIdH[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 }};

PUBLIC int mismatch23I37[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,   -50,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -110,     0,   -70,     0}
 ,{     0,     0,     0,     0,   -30}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -120,     0,   -70,     0}
 ,{     0,     0,     0,     0,   -30}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    20,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    20,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }
,{{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,    70,    70,    70,    70}
 ,{    70,   -40,    70,     0,    70}
 ,{    70,    70,    70,    70,    40}
 }};
PUBLIC int mismatch23IdH[NBPAIRS+1][5][5] =
{{{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,  -570,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,  -860,     0,  -900,     0}
 ,{     0,     0,     0,     0,  -640}
 }
,{{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0,     0,     0,     0,     0}
 ,{     0, -1090,     0,  -900,     0}
 ,{     0,     0,     0,     0,  -640}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -580,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   -60,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -360,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -580,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   -60,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -360,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }
,{{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,   500,   500,   500,   500}
 ,{   500,  -360,   500,  -400,   500}
 ,{   500,   500,   500,   500,  -140}
 }};

/* mismatch_exterior */
PUBLIC int mismatchExt37[NBPAIRS+1][5][5] =
{{ /* NP.. */
  {   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 },
 { /* CG.. */
  {   -50,  -110,   -50,  -140,   -70}
 ,{  -110,  -110,  -110,  -160,  -110}
 ,{   -70,  -150,   -70,  -150,  -100}
 ,{  -110,  -130,  -110,  -140,  -110}
 ,{   -50,  -150,   -50,  -150,   -70}
 },
 { /* GC.. */
  {   -80,  -140,   -80,  -140,  -100}
 ,{  -100,  -150,  -100,  -140,  -100}
 ,{  -110,  -150,  -110,  -150,  -140}
 ,{  -100,  -140,  -100,  -160,  -100}
 ,{   -80,  -150,   -80,  -150,  -120}
 },
 { /* GU.. */
  {   -50,   -80,   -50,   -50,   -50}
 ,{   -50,  -100,   -70,   -50,   -70}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -110,   -70,   -80,   -70}
 ,{   -50,   -80,   -50,   -80,   -50}
 },
 { /* UG.. */
  {   -30,   -30,   -60,   -60,   -60}
 ,{   -30,   -30,   -60,   -60,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -60,  -100,   -70,  -100,   -60}
 },
 { /* AU.. */
  {   -50,   -80,   -50,   -80,   -50}
 ,{   -70,  -100,   -70,  -110,   -70}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -110,   -70,  -120,   -70}
 ,{   -50,   -80,   -50,   -80,   -50}
 },
 { /* UA.. */
  {   -60,   -80,   -60,   -80,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -70,  -100,   -70,  -100,   -80}
 },
 { /* NN.. */
  {   -30,   -30,   -50,   -50,   -50}
 ,{   -30,   -30,   -60,   -50,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -60,   -80,   -60,   -80,   -60}
 ,{   -50,   -80,   -50,   -80,   -50}
 }};

/* mismatch_exterior_enthalpies */
PUBLIC int mismatchExtdH[NBPAIRS+1][5][5] =
{{ /* NP.. */
  {   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 ,{   INF,   INF,   INF,   INF,   INF}
 },
 { /* CG.. */
  {    50,  -400,    50,  -400,   -30}
 ,{  -520,  -520,  -720,  -710,  -720}
 ,{    50,  -400,    50,  -400,   -30}
 ,{  -560,  -560,  -720,  -620,  -720}
 ,{  -400,  -400,  -420,  -400,  -500}
 },
 { /* GC.. */
  {  -270,  -560,  -270,  -560,  -530}
 ,{  -570,  -910,  -570,  -820,  -570}
 ,{  -340,  -560,  -340,  -560,  -530}
 ,{  -560,  -560,  -570,  -920,  -570}
 ,{  -270,  -560,  -270,  -560,  -860}
 },
 { /* GU.. */
  {   310,  -480,  -180,   310,   140}
 ,{   310,  -480,  -430,   310,  -430}
 ,{  -140,  -630,  -510,  -630,  -140}
 ,{  -150,  -890,  -430,  -150,  -430}
 ,{   140,  -630,  -180,  -630,   140}
 },
 { /* UG.. */
  {   600,   200,   600,   200,   460}
 ,{   -60,  -340,  -230,   -60,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -230,  -350,  -230,  -350,  -230}
 ,{   200,   200,   -30,   200,   160}
 },
 { /* AU.. */
  {   140,  -400,  -180,  -380,   140}
 ,{  -380,  -400,  -430,  -380,  -430}
 ,{  -140,  -630,  -510,  -630,  -140}
 ,{  -430,  -890,  -430,  -890,  -430}
 ,{   140,  -630,  -180,  -630,   140}
 },
 { /* UA.. */
  {   600,   200,   600,   200,   460}
 ,{  -230,  -390,  -230,  -310,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -230,  -350,  -230,  -350,  -230}
 ,{   200,   200,   -30,   200,  -170}
 },
 { /* NN.. */
  {   600,   200,   600,   310,   460}
 ,{   310,  -340,  -230,   310,  -230}
 ,{   600,   200,   600,   200,   460}
 ,{  -150,  -350,  -230,  -150,  -230}
 ,{   200,   200,   -30,   200,   160}
 }};

/* dangle5 */
PUBLIC int dangle5_37[NBPAIRS+1][5] =
{ /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {   -10,   -50,   -30,   -20,   -10},
/* GC */ {    -0,   -20,   -30,    -0,    -0},
/* GU */ {   -20,   -30,   -30,   -40,   -20},
/* UG */ {   -10,   -30,   -10,   -20,   -20},
/* AU */ {   -20,   -30,   -30,   -40,   -20},
/* UA */ {   -10,   -30,   -10,   -20,   -20},
/* NN */ {    -0,   -20,   -10,    -0,    -0}
};

/* dangle3 */
PUBLIC int dangle3_37[NBPAIRS+1][5] =
{ /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {   -40,  -110,   -40,  -130,   -60},
/* GC */ {   -80,  -170,   -80,  -170,  -120},
/* GU */ {   -10,   -70,   -10,   -70,   -10},
/* UG */ {   -50,   -80,   -50,   -80,   -60},
/* AU */ {   -10,   -70,   -10,   -70,   -10},
/* UA */ {   -50,   -80,   -50,   -80,   -60},
/* NN */ {   -10,   -70,   -10,   -70,   -10}
};

/* dangle5_enthalpies */
PUBLIC int dangle5_dH[NBPAIRS+1][5] =
{ /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {   330,  -240,   330,    80,  -140},
/* GC */ {    70,  -160,    70,  -460,   -40},
/* GU */ {   310,   160,   220,    70,   310},
/* UG */ {   690,   -50,   690,    60,    60},
/* AU */ {   310,   160,   220,    70,   310},
/* UA */ {   690,   -50,   690,    60,    60},
/* NN */ {   690,   160,   690,    80,   310}
};

/* dangle3_enthalpies */
PUBLIC int dangle3_dH[NBPAIRS+1][5] =
{ /*           N      A      C      G      U */
/* NP */ {   INF,   INF,   INF,   INF,   INF},
/* CG */ {  -280,  -740,  -280,  -640,  -360},
/* GC */ {  -410,  -900,  -410,  -860,  -750},
/* GU */ {   -70,  -570,   -70,  -580,  -220},
/* UG */ {   -90,  -490,   -90,  -550,  -230},
/* AU */ {   -70,  -570,   -70,  -580,  -220},
/* UA */ {   -90,  -490,   -90,  -550,  -230},
/* NN */ {   -70,  -490,   -70,  -550,  -220}
};

PUBLIC char Triloops[241] =
  "CAACG "
  "GUUAC "
;
PUBLIC int Triloop37[40] = {   680,   690};
PUBLIC int TriloopdH[40] = {  2370,  1080};

PUBLIC char Tetraloops[281] =
  "CAACGG "
  "CCAAGG "
  "CCACGG "
  "CCCAGG "
  "CCGAGG "
  "CCGCGG "
  "CCUAGG "
  "CCUCGG "
  "CUAAGG "
  "CUACGG "
  "CUCAGG "
  "CUCCGG "
  "CUGCGG "
  "CUUAGG "
  "CUUCGG "
  "CUUUGG "
;
PUBLIC int Tetraloop37[40] = {   550,   330,   370,   340,   350,   360,   370,   250,   360,   280,   370,   270,   280,   350,   370,   370};
PUBLIC int TetraloopdH[40] = {   690, -1030,  -330,  -890,  -660,  -750,  -350, -1390,  -760, -1070,  -660, -1290, -1070,  -620, -1530,  -680};

PUBLIC char Hexaloops[361] =
  "ACAGUACU "
  "ACAGUGAU "
  "ACAGUGCU "
  "ACAGUGUU "
;
PUBLIC int Hexaloop37[40] = {   280,   360,   290,   180};
PUBLIC int HexaloopdH[40] = { -1680, -1140, -1280, -1540};

#include "intl11.h"
#include "intl11dH.h"
#include "intl21.h"
#include "intl21dH.h"
#include "intl22.h"
#include "intl22dH.h"

