/* unit test for function vrna_eval_structure_pt*/

#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* malloc, free, rand */
#include <math.h>


#include <ViennaRNA/utils.h>
#include <ViennaRNA/exterior_loops.h>
#include <ViennaRNA/params.h>
#include <ViennaRNA/interior_loops.h>
#include <ViennaRNA/hairpin_loops.h>
#include <ViennaRNA/multibranch_loops.h>
#include "ViennaRNA/energy_par.h"
#include <ViennaRNA/data_structures.h>
#include <ViennaRNA/constraints.h>

#include "ViennaRNA/structure_utils.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/loop_energies.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/eval.h"





typedef struct rnaStr
{
    char *sequence;
    char *structure;
    int energy;
    int dangle;
}rnaStr;

#test eval_structure
{
    struct rnaStr strList[15];

    struct rnaStr str;
    int count=0;
{/*definition of structures*/
    /*################################################
     other contributations(loops, stems, ..)
 *  (2,10) GC-CA TM   =                             -110
 *  (3,10)(4,9) stack =                             -330
 *  (4,9) Hairpin = +560 + C-G -AA TM =-150  =      +410
 *  (11,18) same stack                              -330
 *  (12,17) same hairpin                            +410
 *  (2,19) CG-CG TM         =                       -100
 *  (3,10) CG-CC TM         =                       -70
 *  (11,18) CG-GG TM        =                       -140
 *  930 + 3 * -90           =                       +660
 *   energyManuelMultiDangle2 =                     +350

    Overal =                                        +400
 */

    str.sequence=" ACCCAAAAGGCCAAAAGGGC";
    str.structure=".(((....))((....))).";
    str.dangle=2;
    str.energy=400;

    strList[count]=str;
    count++;

    /*################################################

     exterior (24,2) GC-AA TM    -110
     (4,11)CG-AU WC        -210
     (5,10)AU-AA TM        -80
     (5,10) AU closing        50
     hairpin 4            560
     (14,21)AU-CG WC        -220
     CG-AA TM            -150
     hairpin 4            560

     Multiloop(Dangle 2)
    a  + b * branch
    930 + 3 * -90        = 660
    (2,24)CG-AC TM        -150
    (11,4) GC-UA TM        -110
    (21,14) UA-CU TM    -50
    (21,14) AU closure    50
    MUltiloop=        400
     Overall:            800

 */

    str.sequence= "ACACAAAAAUGUUACAAAAGUCCGA";
    str.structure=".(.((....))..((....))..).";
    str.dangle=2;
    str.energy=800;

    strList[count]=str;
    count++;

    /*################################################
     exterior (24,2) UA-AA TM    -100
     (24,2) AU closure         50
     (4,11)CG-AU WC        -210
     (5,10)AU-AA TM        -80
     (5,10) AU closing        50
     hairpin 4            560
     (14,21)AU-CG WC        -220
     CG-AA TM            -150
     hairpin 4            560
                460
     Multiloop(Dangle 2)
    a  + b * branch
    930 + 3 * -90        = 660
    (2,24)AU-AC TM        -100
    (2,24) AU closing    50
    (11,4) GC-UA TM        -110
    (21,14) UA-CU TM    -50
    (21,14) AU closure    50
    MUltiloop=        500
     Overall:            960

 */

    str.sequence= "AAACAAAAAUGUUACAAAAGUCCUA";
    str.structure=".(.((....))..((....))..).";
    str.dangle=2;
    str.energy=960;

    strList[count]=str;
    count++;
    /*################################################
     exterior (24,2) UG-AA TM    -100
     (24,2) GU closure         50
     (4,11)CG-AU WC        -210
     (5,10)AU-AA TM        -80
     (5,10) AU closing        50
     hairpin 4            560
     (14,21)AU-CG WC        -220
     CG-AA TM            -150
     hairpin 4            560
                460
     Multiloop(Dangle 2)
    a  + b * branch
    930 + 3 * -90        = 660
    (2,24)GU-AC TM        -100
    (2,24) GU closing    50
    (11,4) GC-UA TM        -110
    (21,14) UA-CU TM    -50
    (21,14) AU closure    50
    MUltiloop=        500
     Overall:            960

 */

    str.sequence= "AGACAAAAAUGUUACAAAAGUCCUA";
    str.structure=".(.((....))..((....))..).";
    str.dangle=2;
    str.energy=960;

    strList[count]=str;
    count++;
    /*################################################
     other contributations(loops, stems, ..)
 *  (2,10) CG-AC TM   =                             -110
 *  (3,10)(4,9) stack =                             -330
 *  (4,9) Hairpin = +560 + C-G -AA TM =-150  =      +410
 *  (11,18) same stack                              -330
 *  (12,17) same hairpin                            +410
 *                                                  +50
 *  Multiloop Dangle option 3
 *   (3,10),(11,18) GC-CG WC pair(flush coaxial) = -340
 *   930 + 3 * -90           =                     +660
 *   energyManuelMultiDangle3                    = +320
    Overall=                                      +370*/

    str.sequence="ACCCAAAAGGCCAAAAGGGC";
    str.structure=".(((....))((....))).";
    str.dangle=3;
    str.energy=370;

    strList[count]=str;
    count++;
    /*################################################

     exterior  GC-AA TM        -110
     1.Branch
     (4,11)CG-AU WC        -210
     (5,10)AU-AA TM        -80
     (5,10) AU closing        50
     hairpin 4            560
     2.Branch
     (14,21)AU-CG WC        -220
     CG-AA TM            -150
     hairpin 4            560
     3.Branch
     CG-CG WC            -330
     CG-AA TM            -150
     hairpin 4            560
                480
     Multiloop(Dangle 2)
    a  + b * branch
    930 + 4 * -90        = 570
    (2,24)CG-AC TM        -150
    (11,4) GC-UA TM        -110
    (21,14) UA-CU TM    -50
    (21,14) AU closure    50
    GC-CU TM(at 3.branch)    -50
    MUltiloop=        260
     Overall:            740

 */

    str.sequence= "ACACAAAAAUGUUACAAAAGUCCAAAAGGCGA";
    str.structure=".(.((....))..((....))((....)).).";
    str.dangle=2;
    str.energy=740;

    strList[count]=str;
    count++;
    /*################################################*/



    str.sequence=" CAAAAAAG";
    str.structure="(......)";
    str.dangle=2;
    str.energy=-150+540;

    strList[count]=str;
    count++;

    /*################################################*/
    str.sequence=" CCAAAAAAGG";
    str.structure="((......))";
    str.dangle=2;
    str.energy=-330-150+540;

    strList[count]=str;
    count++;

    /*################################################*/

    str.sequence=" CAAAAAAAUG";
    str.structure="((......))";
    str.dangle=2;
    str.energy=-210-80+50+540;

    strList[count]=str;
    count++;

    /*################################################*/
    str.sequence=" GAAAAAACCC";
    str.structure="(........)";
    str.dangle=2;
    str.energy=550-150;

    strList[count]=str;
    count++;

    /*################################################*/
    str.sequence=" CACAAAAAAGAG";
    str.structure="(.(......).)";
    str.dangle=2;
    str.energy=90-150+540;

    strList[count]=str;
    count++;
    /*################################################*/


    str.sequence=" CAGCAAAAAAGAUG";
    str.structure="((.(......).))";
    str.dangle=2;
    str.energy=-210+120-150+540;

    strList[count]=str;
    count++;
    /*################################################*/

    str.sequence=" AGCAAAAAAGAU";
    str.structure="(.(......).)";
    str.dangle=2;
    str.energy=50+120-150+540;

    strList[count]=str;
    count++;
    /*################################################*/



    str.sequence=" CUAAAAAAAUCAG";
    str.structure="(.(......)..)";
    str.dangle=2;
    str.energy=300-80+50+540;

    strList[count]=str;
    count++;
    /*################################################*/
    str.sequence=" UGCAAAAAAAUGA";
    str.structure="(..(......).)";
    str.dangle=2;
    str.energy=50+260-80+50+540;

    strList[count]=str;
    count++;

     /*################################################*/

}


    char s[100]="";
    vrna_fold_compound_t *vc;
    short * pairtable;
    int i;
    for(i=0;i< sizeof(strList)/sizeof(strList[0]);i++)
    {
        vc = vrna_fold_compound(strList[i].sequence, NULL, VRNA_OPTION_MFE | VRNA_OPTION_EVAL_ONLY);
        vc->params->model_details.dangles=strList[i].dangle;
        pairtable = vrna_ptable(strList[i].structure);
        int vrnaEnergy=vrna_eval_structure_pt(vc,pairtable);

        ck_assert_msg(strList[i].energy==vrnaEnergy,
                      "\n structure: %s   sequence: %s   manually = %i , vRNA =  %i\n",
                      strList[i].structure, strList[i].sequence, strList[i].energy, vrnaEnergy);
        free(pairtable);
        vrna_fold_compound_free(vc);
    }

}
