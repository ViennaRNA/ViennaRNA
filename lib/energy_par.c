

/*
    Automatically generated using the TurnerParser
    TurnerParser (c) 2008,2009
      Christian Hoener zu Siederdissen, TBI Vienna
      choener (at) tbi.univie.ac.at
*/

/*
     Current free energy parameters are summarized in:

     D.H.Mathews, J. Sabina, M. ZUker, D.H. Turner
     "Expanded sequence dependence of thermodynamic parameters improves
     prediction of RNA secondary structure"
     JMB, 288, pp 911-940, 1999

     Enthalpies taken from:

     A. Walter, D Turner, J Kim, M Lyttle, P M"uller, D Mathews, M Zuker
     "Coaxial stckaing of helices enhances binding of oligoribonucleotides.."
     PNAS, 91, pp 9218-9222, 1994

     D.H. Turner, N. Sugimoto, and S.M. Freier.
     "RNA Structure Prediction",
     Ann. Rev. Biophys. Biophys. Chem. 17, 167-192, 1988.

     John A.Jaeger, Douglas H.Turner, and Michael Zuker.
     "Improved predictions of secondary structures for RNA",
     PNAS, 86, 7706-7710, October 1989.

     L. He, R. Kierzek, J. SantaLucia, A.E. Walter, D.H. Turner
     "Nearest-Neughbor Parameters for GU Mismatches...."
     Biochemistry 1991, 30 11124-11132

     A.E. Peritz, R. Kierzek, N, Sugimoto, D.H. Turner
     "Thermodynamic Study of Internal Loops in Oligoribonucleotides..."
     Biochemistry 1991, 30, 6428--6435
*/



#include "energy_const.h"
/*@unused@*/
static char rcsid[] = "$Id: energy_par.c,v 1.6 2004/08/12 12:11:57 ivo Exp $";

#define NST 0     /* Energy for nonstandard stacked pairs */
#define DEF -50   /* Default terminal mismatch, used for I */
                  /* and any non_pairing bases */
#define NSM 0     /* terminal mismatch for non standard pairs */

#define PUBLIC

PUBLIC double Tmeasure = 37+K0;  /* temperature of param measurements */



//PUBLIC double lxc37=1.079;
PUBLIC double lxc37=107.856;

/* Multi loops */
PUBLIC int ML_BASE37=0;
PUBLIC int ML_closing37=930;
PUBLIC int ML_intern37=-90;
PUBLIC int ML_BASEdH=0;
PUBLIC int ML_closingdH=3000;
PUBLIC int ML_interndH=-220;

/* Ninio */
PUBLIC int MAX_NINIO=300;
PUBLIC int ninio37=60;
PUBLIC int niniodH=320;

/* terminal AU */
PUBLIC int TerminalAU37=50;
PUBLIC int TerminalAUdH=370;

/* duplex */
PUBLIC int DuplexInit37=410;
PUBLIC int DuplexInitdH=360;

PUBLIC int stack37[NBPAIRS+1][NBPAIRS+1] = {
{    INF,    INF,    INF,    INF,    INF,    INF,    INF,    INF},
{    INF,   -240,   -330,   -210,   -140,   -210,   -210,    NST},
{    INF,   -330,   -340,   -250,   -150,   -220,   -240,    NST},
{    INF,   -210,   -250,    130,    -50,   -140,   -130,    NST},
{    INF,   -140,   -150,    -50,     30,    -60,   -100,    NST},
{    INF,   -210,   -220,   -140,    -60,   -110,    -90,    NST},
{    INF,   -210,   -240,   -130,   -100,    -90,   -130,    NST},
{    INF,    NST,    NST,    NST,    NST,    NST,    NST,    NST}};

PUBLIC int stackdH[NBPAIRS+1][NBPAIRS+1] = {
{    INF,    INF,    INF,    INF,    INF,    INF,    INF,    INF},
{    INF,  -1060,  -1340,  -1210,   -560,  -1050,  -1040,    NST},
{    INF,  -1340,  -1490,  -1260,   -830,  -1140,  -1240,    NST},
{    INF,  -1210,  -1260,  -1460,  -1350,   -880,  -1280,    NST},
{    INF,   -560,   -830,  -1350,   -930,   -320,   -700,    NST},
{    INF,  -1050,  -1140,   -880,   -320,   -940,   -680,    NST},
{    INF,  -1040,  -1240,  -1280,   -700,   -680,   -770,    NST},
{    INF,    NST,    NST,    NST,    NST,    NST,    NST,    NST}};

/* loops of different kinds */
PUBLIC int hairpin37[31] = {    INF,    INF,    INF,    540,    560,    570,    540,    600,    550,    640,    650,    660,    670,    680,    690,    690,    700,    710,    710,    720,    720,    730,    730,    740,    740,    750,    750,    750,    760,    760,    770};
PUBLIC int hairpindH[31] = {    INF,    INF,    INF,    130,    480,    360,   -290,    130,   -290,    500,    500,    500,    500,    500,    500,    500,    500,    500,    500,    500,    500,    500,    500,    500,    500,    500,    500,    500,    500,    500,    500};
PUBLIC int bulge37[31] = {    INF,    380,    280,    320,    360,    400,    440,    460,    470,    480,    490,    500,    510,    520,    530,    540,    540,    550,    550,    560,    570,    570,    580,    580,    580,    590,    590,    600,    600,    600,    610};
PUBLIC int bulgedH[31] = {    INF,   1060,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710,    710};
PUBLIC int internal_loop37[31] = {    INF,    INF,    INF,    INF,    110,    200,    200,    210,    230,    240,    250,    260,    270,    280,    290,    290,    300,    310,    310,    320,    330,    330,    340,    340,    350,    350,    350,    360,    360,    370,    370};
PUBLIC int internal_loopdH[31] = {    INF,    INF,    INF,    INF,   -720,   -680,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130,   -130};

PUBLIC int mismatchI37[NBPAIRS+1][5][5] = 
{{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,      0,      0,    -80,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,   -100,      0,   -100,      0}
 ,{    INF,      0,      0,      0,    -60}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,      0,      0,    -80,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,   -100,      0,   -100,      0}
 ,{    INF,      0,      0,      0,    -60}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,     70,     70,    -10,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,    -30,     70,    -30,     70}
 ,{    INF,     70,     70,     70,     10}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,     70,     70,    -10,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,    -30,     70,    -30,     70}
 ,{    INF,     70,     70,     70,     10}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,     70,     70,    -10,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,    -30,     70,    -30,     70}
 ,{    INF,     70,     70,     70,     10}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,     70,     70,    -10,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,    -30,     70,    -30,     70}
 ,{    INF,     70,     70,     70,     10}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
                                          }
                                          }
;

PUBLIC int mismatchIdH[NBPAIRS+1][5][5] = 
{{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,      0,      0,   -340,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,   -760,      0,    280,      0}
 ,{    INF,      0,      0,      0,   -580}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,      0,      0,   -340,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,   -760,      0,    280,      0}
 ,{    INF,      0,      0,      0,   -580}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    500,    500,    170,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,   -260,    500,    790,    500}
 ,{    INF,    500,    500,    500,    -80}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    500,    500,    170,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,   -260,    500,    790,    500}
 ,{    INF,    500,    500,    500,    -80}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    500,    500,    170,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,   -260,    500,    790,    500}
 ,{    INF,    500,    500,    500,    -80}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    500,    500,    170,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,   -260,    500,    790,    500}
 ,{    INF,    500,    500,    500,    -80}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
                                          }
                                          }
;

PUBLIC int mismatchH37[NBPAIRS+1][5][5] = 
{{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -150,   -150,   -140,   -150}
 ,{    INF,   -100,   -110,   -100,    -80}
 ,{    INF,   -230,   -150,   -240,   -150}
 ,{    INF,   -100,   -140,   -100,   -210}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -110,   -150,   -130,   -150}
 ,{    INF,   -110,    -70,   -110,    -50}
 ,{    INF,   -250,   -150,   -220,   -150}
 ,{    INF,   -110,   -100,   -110,   -160}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,     20,    -50,    -30,    -50}
 ,{    INF,    -10,    -20,    -10,    -20}
 ,{    INF,   -100,    -50,   -110,    -50}
 ,{    INF,    -10,    -30,    -10,   -100}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    -50,    -30,    -60,    -30}
 ,{    INF,    -20,    -10,    -20,      0}
 ,{    INF,    -90,    -30,   -110,    -30}
 ,{    INF,    -20,    -10,    -20,    -90}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    -30,    -50,    -30,    -50}
 ,{    INF,    -10,    -20,    -10,    -20}
 ,{    INF,   -120,    -50,   -110,    -50}
 ,{    INF,    -10,    -30,    -10,   -120}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    -50,    -30,    -50,    -30}
 ,{    INF,    -20,    -10,    -20,      0}
 ,{    INF,   -150,    -30,   -150,    -30}
 ,{    INF,    -20,    -10,    -20,    -90}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
                                          }
                                          }
;

PUBLIC int mismatchHdH[NBPAIRS+1][5][5] = 
{{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -910,   -560,   -560,   -560}
 ,{    INF,   -570,   -340,   -570,   -270}
 ,{    INF,  -1400,    560,   -920,   -560}
 ,{    INF,   -570,   -530,   -570,  -1440}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -520,   -400,   -560,   -400}
 ,{    INF,   -720,     50,   -720,   -420}
 ,{    INF,  -1290,   -400,   -620,   -400}
 ,{    INF,   -720,    -30,   -720,  -1080}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,     30,    570,     20,    570}
 ,{    INF,    140,    970,    140,    340}
 ,{    INF,   -270,    570,     20,    570}
 ,{    INF,    140,    830,    140,    -50}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -110,   -260,   -520,   -260}
 ,{    INF,    -60,   -140,    -60,    190}
 ,{    INF,    100,   -260,    220,   -260}
 ,{    INF,    -60,    230,    -60,    -70}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    -20,    570,     20,    570}
 ,{    INF,    140,    970,    140,    340}
 ,{    INF,   -520,    570,     20,    570}
 ,{    INF,    140,    830,    140,   -380}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    -30,   -260,   -520,   -260}
 ,{    INF,    -60,   -140,    -60,    190}
 ,{    INF,   -590,   -260,   -520,   -260}
 ,{    INF,    -60,    230,    -60,    -70}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
                                          }
                                          }
;

PUBLIC int mismatchM37[NBPAIRS+1][5][5] = 
{{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -150,   -150,   -140,   -150}
 ,{    INF,   -100,   -110,   -100,    -80}
 ,{    INF,   -140,   -150,   -160,   -150}
 ,{    INF,   -100,   -140,   -100,   -120}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -110,   -150,   -130,   -150}
 ,{    INF,   -110,    -70,   -110,    -50}
 ,{    INF,   -160,   -150,   -140,   -150}
 ,{    INF,   -110,   -100,   -110,    -70}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    -30,   -100,    -80,   -100}
 ,{    INF,    -60,    -70,    -60,    -70}
 ,{    INF,    -60,   -100,    -80,   -100}
 ,{    INF,    -60,    -80,    -60,    -60}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -100,    -80,   -110,    -80}
 ,{    INF,    -70,    -60,    -70,    -50}
 ,{    INF,    -50,    -80,    -80,    -80}
 ,{    INF,    -70,    -60,    -70,    -50}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    -80,   -100,    -80,   -100}
 ,{    INF,    -60,    -70,    -60,    -70}
 ,{    INF,    -80,   -100,    -80,   -100}
 ,{    INF,    -60,    -80,    -60,    -80}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -100,    -80,   -110,    -80}
 ,{    INF,    -70,    -60,    -70,    -50}
 ,{    INF,   -110,    -80,   -120,    -80}
 ,{    INF,    -70,    -60,    -70,    -50}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
                                          }
                                          }
;

PUBLIC int mismatchMdH[NBPAIRS+1][5][5] = 
{{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -910,   -560,   -560,   -560}
 ,{    INF,   -570,   -340,   -570,   -270}
 ,{    INF,   -820,   -560,   -920,   -560}
 ,{    INF,   -570,   -530,   -570,   -860}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -520,   -400,   -560,   -400}
 ,{    INF,   -720,     50,   -720,   -420}
 ,{    INF,   -710,   -400,   -620,   -400}
 ,{    INF,   -720,    -30,   -720,   -500}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -340,    200,   -350,    200}
 ,{    INF,   -230,    600,   -230,    -30}
 ,{    INF,    -60,    200,   -350,    200}
 ,{    INF,   -230,    460,   -230,    160}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -480,   -630,   -890,   -630}
 ,{    INF,   -430,   -510,   -430,   -180}
 ,{    INF,    310,   -630,   -150,   -630}
 ,{    INF,   -430,   -140,   -430,    140}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -390,    200,   -350,    200}
 ,{    INF,   -230,    600,   -230,    -30}
 ,{    INF,   -310,    200,   -350,    200}
 ,{    INF,   -230,    460,   -230,   -170}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -400,   -630,   -890,   -630}
 ,{    INF,   -430,   -510,   -430,   -180}
 ,{    INF,   -380,   -630,   -890,   -630}
 ,{    INF,   -430,   -140,   -430,    140}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
                                          }
                                          }
;

PUBLIC int mismatch1nI37[NBPAIRS+1][5][5] = 
{{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,      0,      0,      0,      0}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,      0,      0,      0,      0}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,     70,     70,     70,     70}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,     70,     70,     70,     70}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,     70,     70,     70,     70}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,     70,     70,     70,     70}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
                                          }
                                          }
;

PUBLIC int mismatch1nIdH[NBPAIRS+1][5][5] = 
{{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,      0,      0,      0,      0}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,      0,      0,      0,      0}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,    500,    500,    500,    500}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,    500,    500,    500,    500}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,    500,    500,    500,    500}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,    500,    500,    500,    500}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
                                          }
                                          }
;

PUBLIC int mismatch23I37[NBPAIRS+1][5][5] = 
{{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,      0,      0,    -50,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,   -110,      0,    -70,      0}
 ,{    INF,      0,      0,      0,    -30}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,   -120,      0,    -70,      0}
 ,{    INF,      0,      0,      0,    -30}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,    -40,     70,      0,     70}
 ,{    INF,     70,     70,     70,     40}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,     70,     70,     20,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,    -40,     70,      0,     70}
 ,{    INF,     70,     70,     70,     40}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,    -40,     70,      0,     70}
 ,{    INF,     70,     70,     70,     40}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,     70,     70,     20,     70}
 ,{    INF,     70,     70,     70,     70}
 ,{    INF,    -40,     70,      0,     70}
 ,{    INF,     70,     70,     70,     40}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
                                          }
                                          }
;

PUBLIC int mismatch23IdH[NBPAIRS+1][5][5] = 
{{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,      0,      0,   -570,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,   -860,      0,   -900,      0}
 ,{    INF,      0,      0,      0,   -640}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,      0,      0,      0,      0}
 ,{    INF,  -1090,      0,   -900,      0}
 ,{    INF,      0,      0,      0,   -640}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,   -580,    500,   -400,    500}
 ,{    INF,    500,    500,    500,   -140}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    500,    500,    -60,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,   -360,    500,   -400,    500}
 ,{    INF,    500,    500,    500,   -140}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,   -580,    500,   -400,    500}
 ,{    INF,    500,    500,    500,   -140}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    500,    500,    -60,    500}
 ,{    INF,    500,    500,    500,    500}
 ,{    INF,   -360,    500,   -400,    500}
 ,{    INF,    500,    500,    500,   -140}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
                                          }
                                          }
;

PUBLIC int mismatchExt37[NBPAIRS+1][5][5] = 
{{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -150,   -150,   -140,   -150}
 ,{    INF,   -100,   -110,   -100,    -80}
 ,{    INF,   -140,   -150,   -160,   -150}
 ,{    INF,   -100,   -140,   -100,   -120}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -110,   -150,   -130,   -150}
 ,{    INF,   -110,    -70,   -110,    -50}
 ,{    INF,   -160,   -150,   -140,   -150}
 ,{    INF,   -110,   -100,   -110,    -70}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    -30,   -100,    -80,   -100}
 ,{    INF,    -60,    -70,    -60,    -70}
 ,{    INF,    -60,   -100,    -80,   -100}
 ,{    INF,    -60,    -80,    -60,    -60}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -100,    -80,   -110,    -80}
 ,{    INF,    -70,    -60,    -70,    -50}
 ,{    INF,    -50,    -80,    -80,    -80}
 ,{    INF,    -70,    -60,    -70,    -50}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    -80,   -100,    -80,   -100}
 ,{    INF,    -60,    -70,    -60,    -70}
 ,{    INF,    -80,   -100,    -80,   -100}
 ,{    INF,    -60,    -80,    -60,    -80}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -100,    -80,   -110,    -80}
 ,{    INF,    -70,    -60,    -70,    -50}
 ,{    INF,   -110,    -80,   -120,    -80}
 ,{    INF,    -70,    -60,    -70,    -50}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
                                          }
                                          }
;

PUBLIC int mismatchExtdH[NBPAIRS+1][5][5] = 
{{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    INF,    INF,    INF,    INF}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -910,   -560,   -560,   -560}
 ,{    INF,   -570,   -340,   -570,   -270}
 ,{    INF,   -820,   -560,   -920,   -560}
 ,{    INF,   -570,   -530,   -570,   -860}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -520,   -400,   -560,   -400}
 ,{    INF,   -720,     50,   -720,   -420}
 ,{    INF,   -710,   -400,   -620,   -400}
 ,{    INF,   -720,    -30,   -720,   -500}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -340,    200,   -350,    200}
 ,{    INF,   -230,    600,   -230,    -30}
 ,{    INF,    -60,    200,   -350,    200}
 ,{    INF,   -230,    460,   -230,    160}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -480,   -630,   -890,   -630}
 ,{    INF,   -430,   -510,   -430,   -180}
 ,{    INF,    310,   -630,   -150,   -630}
 ,{    INF,   -430,   -140,   -430,    140}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -390,    200,   -350,    200}
 ,{    INF,   -230,    600,   -230,    -30}
 ,{    INF,   -310,    200,   -350,    200}
 ,{    INF,   -230,    460,   -230,   -170}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,   -400,   -630,   -890,   -630}
 ,{    INF,   -430,   -510,   -430,   -180}
 ,{    INF,   -380,   -630,   -890,   -630}
 ,{    INF,   -430,   -140,   -430,    140}
                                          }
,{{    INF,    INF,    INF,    INF,    INF}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
 ,{    INF,    NST,    NST,    NST,    NST}
                                          }
                                          }
;

PUBLIC int dangle3_37[NBPAIRS+1][5] = 
{{    INF,    INF,    INF,    INF,    INF}
,{    INF,   -170,    -80,   -170,   -120}
,{    INF,   -110,    -40,   -130,    -60}
,{    INF,    -80,    -50,    -80,    -60}
,{    INF,    -70,    -10,    -70,    -10}
,{    INF,    -80,    -50,    -80,    -60}
,{    INF,    -70,    -10,    -70,    -10}
,{    INF,    NST,    NST,    NST,    NST}
                                         }
;

PUBLIC int dangle3_dH[NBPAIRS+1][5] = 
{{    INF,    INF,    INF,    INF,    INF}
,{    INF,   -900,   -410,   -860,   -750}
,{    INF,   -740,   -280,   -640,   -360}
,{    INF,   -490,    -90,   -550,   -230}
,{    INF,   -570,    -70,   -580,   -220}
,{    INF,   -490,    -90,   -550,   -230}
,{    INF,   -570,    -70,   -580,   -220}
,{    INF,    NST,    NST,    NST,    NST}
                                         }
;

PUBLIC int dangle5_37[NBPAIRS+1][5] = 
{{    INF,    INF,    INF,    INF,    INF}
,{    INF,    -20,    -30,      0,      0}
,{    INF,    -50,    -30,    -20,    -10}
,{    INF,    -30,    -10,    -20,    -20}
,{    INF,    -30,    -30,    -40,    -20}
,{    INF,    -30,    -10,    -20,    -20}
,{    INF,    -30,    -30,    -40,    -20}
,{    INF,    NST,    NST,    NST,    NST}
                                         }
;

PUBLIC int dangle5_dH[NBPAIRS+1][5] = 
{{    INF,    INF,    INF,    INF,    INF}
,{    INF,   -900,   -410,   -860,   -750}
,{    INF,   -740,   -280,   -640,   -360}
,{    INF,   -490,    -90,   -550,   -230}
,{    INF,   -570,    -70,   -580,   -220}
,{    INF,   -490,    -90,   -550,   -230}
,{    INF,   -570,    -70,   -580,   -220}
,{    INF,    NST,    NST,    NST,    NST}
                                         }
;

PUBLIC char Triloops[241] =
  "CAACG "
  "GUUAC "
;
PUBLIC int Triloop37[40] = {    680,    690};
PUBLIC int TriloopdH[40] = {   2370,   1080};


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
PUBLIC int Tetraloop37[40] = {    550,    330,    370,    340,    350,    360,    370,    250,    360,    280,    370,    270,    280,    350,    370,    370};
PUBLIC int TetraloopdH[40] = {    690,  -1030,   -330,   -890,   -660,   -750,   -350,  -1390,   -760,  -1070,   -660,  -1290,  -1070,   -620,  -1530,   -680};


PUBLIC char Hexaloops[361] =
  "ACAGUACU "
  "ACAGUGAU "
  "ACAGUGCU "
  "ACAGUGUU "
;
PUBLIC int Hexaloop37[40] = {    280,    360,    290,    180};
PUBLIC int HexaloopdH[40] = {  -1680,  -1140,  -1280,  -1540};


#include "intl11.h"
#include "intl11dH.h"
#include "intl21.h"
#include "intl21dH.h"
#include "intl22.h"
#include "intl22dH.h"
