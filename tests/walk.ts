
#include <stdio.h>
#include <stdlib.h>
#include <ViennaRNA/walk.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/data_structures.h>

#suite Walks

#test Walk_Gradient

char *sequence =                            "GGGAAACCCAACCUUU";
char *structure =                           ".(.....)........";
char *expectedGradientBasinRepresentative = "(((...))).......";
int expectedLength = 2;
move expectedMoves[2] ={{1,9},{3,7}};
int shifts = 0;
vrna_md_t md;
vrna_md_set_default(&md);
vrna_fold_compound_t *vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
short *pt = vrna_ptable(structure);

//WalkTest
short *basinRep = vrna_ptable_copy(pt);
move *resultMoves = vrna_walk_gradient_path(vc, basinRep, shifts);
char *resultStructure = vrna_db_from_ptable(basinRep);

//check if the final structure is correct.
ck_assert_str_eq( resultStructure, expectedGradientBasinRepresentative);

//check if all moves are in the set.
int movesFound = 0;
for(int i = 0; i < expectedLength;i++){
  for(move *m = resultMoves; m->bpLeft != 0; m++){
    if(expectedMoves[i].bpLeft == m->bpLeft && expectedMoves[i].bpRight == m->bpRight){
      movesFound++;
      break;
    }
  }
}

ck_assert_int_eq(movesFound, expectedLength);


vrna_fold_compound_free(vc);
free(pt);
free(basinRep);
free(resultMoves);
free(resultStructure);

