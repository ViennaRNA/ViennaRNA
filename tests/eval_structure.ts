/* unit test for function structure energy evaluations */

#include <stdio.h>
#include <stdlib.h>


#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/utils/structures.h>
#include "ViennaRNA/eval.h"

typedef struct {
  char  *sequence;
  char  *structure;
  int   energy;
  int   dangle;
} rnaStr;

#test eval_structure
{
  rnaStr  strList[20];

  rnaStr  str;
  int     count = 0;

  {
    str.sequence  = "ACCCAAAAGGCCAAAAGGGC";
    str.structure = ".(((....))((....))).";
    str.dangle    = 2;
    str.energy    = 400;

    strList[count++] = str;

    str.sequence  = "ACACAAAAAUGUUACAAAAGUCCGA";
    str.structure = ".(.((....))..((....))..).";
    str.dangle    = 2;
    str.energy    = 800;

    strList[count++] = str;

    str.sequence  = "AAACAAAAAUGUUACAAAAGUCCUA";
    str.structure = ".(.((....))..((....))..).";
    str.dangle    = 2;
    str.energy    = 960;

    strList[count++] = str;

    str.sequence  = "AGACAAAAAUGUUACAAAAGUCCUA";
    str.structure = ".(.((....))..((....))..).";
    str.dangle    = 2;
    str.energy    = 960;

    strList[count++] = str;

    str.sequence  = "ACCCAAAAGGCCAAAAGGGC";
    str.structure = ".(((....))((....))).";
    str.dangle    = 3;
    str.energy    = 370;

    strList[count++] = str;

    str.sequence  = "ACACAAAAAUGUUACAAAAGUCCAAAAGGCGA";
    str.structure = ".(.((....))..((....))((....)).).";
    str.dangle    = 2;
    str.energy    = 740;

    strList[count++] = str;

    str.sequence  = "CAAAAAAG";
    str.structure = "(......)";
    str.dangle    = 2;
    str.energy    = -150 + 540;

    strList[count++] = str;

    str.sequence  = "CCAAAAAAGG";
    str.structure = "((......))";
    str.dangle    = 2;
    str.energy    = -330 - 150 + 540;

    strList[count++] = str;

    str.sequence  = "CAAAAAAAUG";
    str.structure = "((......))";
    str.dangle    = 2;
    str.energy    = -210 - 80 + 50 + 540;

    strList[count++] = str;

    str.sequence  = "GAAAAAACCC";
    str.structure = "(........)";
    str.dangle    = 2;
    str.energy    = 550 - 150;

    strList[count++] = str;

    str.sequence  = "CACAAAAAAGAG";
    str.structure = "(.(......).)";
    str.dangle    = 2;
    str.energy    = 90 - 150 + 540;

    strList[count++] = str;

    str.sequence  = "CAGCAAAAAAGAUG";
    str.structure = "((.(......).))";
    str.dangle    = 2;
    str.energy    = -210 + 120 - 150 + 540;

    strList[count++] = str;

    str.sequence  = "AGCAAAAAAGAU";
    str.structure = "(.(......).)";
    str.dangle    = 2;
    str.energy    = 50 + 120 - 150 + 540;

    strList[count++] = str;

    str.sequence  = "CUAAAAAAAUCAG";
    str.structure = "(.(......)..)";
    str.dangle    = 2;
    str.energy    = 300 - 80 + 50 + 540;

    strList[count++] = str;

    str.sequence  = "UGCAAAAAAAUGA";
    str.structure = "(..(......).)";
    str.dangle    = 2;
    str.energy    = 50 + 260 - 80 + 50 + 540;

    strList[count++] = str;
  }

  char                  s[100] = "";
  vrna_fold_compound_t  *vc;
  short                 *pairtable;
  int                   i;
  for (i = 0; i < count; i++) {
    vc = vrna_fold_compound(strList[i].sequence,
                            NULL,
                            VRNA_OPTION_DEFAULT);
    vc->params->model_details.dangles = strList[i].dangle;
    pairtable                         = vrna_ptable(strList[i].structure);
    int vrnaEnergy = vrna_eval_structure_pt(vc, pairtable);

    ck_assert_msg(strList[i].energy == vrnaEnergy,
                  "\n structure: %s   sequence: %s   manually = %i , vRNA =  %i\n",
                  strList[i].structure, strList[i].sequence, strList[i].energy, vrnaEnergy);
    free(pairtable);
    vrna_fold_compound_free(vc);
  }
}
