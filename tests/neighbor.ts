#include <ViennaRNA/neighbor.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/data_structures.h>

#suite Neighbor

#test test_vrna_build_neighbors

char *sequence = "GGGAAACCCAACCUUU";
char *structure = ".(.....)........";
int expectedLength = 8;
move expectedNeighbors[8] =
    { { 1, 9 }, { 3, 7 }, { -2, -8 }, { 10, 16 }, { 10, 15 }, { 10, 14 }, { 11, 16 }, { 11, 15 } };
vrna_md_t md;
vrna_md_set_default (&md);
vrna_fold_compound_t *vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
short *pt = vrna_ptable(structure);

int *loopIndices = vrna_loopidx_from_ptable(pt);
int shifts = 0;
move * neighbors = vrna_build_neighbors(vc, pt, shifts);

//check if all neighbors are in the set.
int neighborsFound = 0;
for(int i = 0; i < expectedLength;i++){
  for(move *m = neighbors; m->bpLeft != 0; m++){
    if(expectedNeighbors[i].bpLeft == m->bpLeft && expectedNeighbors[i].bpRight == m->bpRight){
      neighborsFound++;
      break;
    }
  }
}

ck_assert_int_eq(neighborsFound, expectedLength);

vrna_fold_compound_free (vc);
free (pt);
free (loopIndices);
free (neighbors);


#test test_vrna_build_neighbors_with_shifts

char *sequence =  "GAAACCAACCU";
char *structure = "(....).....";
//                "(....)(...)" after move
//                "123456789..
int expectedLength = 6;
int shifts = 1;
move expectedNeighbors[6] =
    { 
        { -1, -6 }, { 7, 11 }, 
        { 1, -5 }, { 1, -9 }, { 1, -10 }, { 1, -11 } /* shifts */
    };
vrna_md_t md;
vrna_md_set_default (&md);
vrna_fold_compound_t *vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
short *pt = vrna_ptable(structure);

int *loopIndices = vrna_loopidx_from_ptable(pt);

move * neighbors = vrna_build_neighbors(vc, pt,shifts);

//check if all neighbors are in the set.
int neighborsFound = 0;
for(int i = 0; i < expectedLength;i++){
  for(move *m = neighbors; m->bpLeft != 0; m++){
    if(expectedNeighbors[i].bpLeft == m->bpLeft && expectedNeighbors[i].bpRight == m->bpRight){
      neighborsFound++;
      break;
    }
  }
}

ck_assert_int_eq(neighborsFound, expectedLength);

vrna_fold_compound_free (vc);
free (pt);
free (loopIndices);
free (neighbors);


#test test_vrna_perform_move

char *sequence = "GGGAAACCCAACCUUU";
char *structure = "................";
short *pt = vrna_ptable(structure);
move m = {2, 8};
vrna_perform_move(&m, pt);
char * resultInsertion = vrna_db_from_ptable(pt);
char *expectedResult = ".(.....)........";
//check if the move was inserted.
ck_assert_str_eq( resultInsertion, expectedResult );

//check if deletion is possible
move m2 = {-2, -8};
vrna_perform_move(&m2, pt);
char * resultDeletion = vrna_db_from_ptable(pt);
ck_assert_str_eq( resultDeletion, structure );

free (resultInsertion);
free (resultDeletion);


#test test_update_loop_indices_deletion

char *sequence = "GGGAAACCCAACCUUU";
char *structure = ".(.....)........";
//test deletion
move m = { -2, -8 };
int expectedLoopIndices[17] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
short *pt = vrna_ptable(structure);
int *loopIndices = vrna_loopidx_from_ptable(pt);

vrna_update_loop_indices(strlen(structure), loopIndices, &m, pt);

int notEqual = 0;
for(int i =0; i < pt[0]; i++){
  if(loopIndices[i] != expectedLoopIndices[i]){
    notEqual = 1;
  }
}
ck_assert_int_eq(notEqual, 0);

free (loopIndices);
free (pt);


#test test_update_loop_indices_insertion
//test insertion
//                "GGGAAACCCAACCUUU";
char *structure = "................";
int expectedLoopIndices1[17] = { 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 };
short *pt = vrna_ptable(structure);
int *loopIndices = vrna_loopidx_from_ptable(pt);
move m1 ={2,8};
vrna_update_loop_indices(strlen(structure), loopIndices, &m1, pt);
int notEqual = 0;
for(int i =0; i < pt[0]; i++){
  if(loopIndices[i] != expectedLoopIndices1[i]){
    notEqual = 1;
  }
}
ck_assert_int_eq(notEqual, 0);

free (loopIndices);
free (pt);

#test test_update_loop_indices_insertion_after_insertion
//test insertion after insertion
//          "GGGAAACCCAACCUUU";
char *structure = ".(.....)........";
//          ".(.....)..(...)." after move
int expectedLoopIndices2[17] = { 2, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 2, 2, 2, 2, 2, 0 };
short *pt = vrna_ptable(structure);
int *loopIndices = vrna_loopidx_from_ptable(pt);
move m2 ={11,15};
vrna_update_loop_indices(strlen(structure), loopIndices, &m2, pt);
int notEqual = 0;
for(int i =0; i < pt[0]; i++){
  if(loopIndices[i] != expectedLoopIndices2[i]){
    notEqual = 1;
  }
}
ck_assert_int_eq(notEqual, 0);

free (loopIndices);
free (pt);


#test test_update_loop_indices_insertion_within_insertion
//test insertion within insertion
//          "GGGAAACCCAACCUUU";
char * structure = ".(.....)........";
//          ".((...))........" after move
int expectedLoopIndices3[17] = { 2, 0, 1, 2, 2, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0 };
short *pt = vrna_ptable(structure);
int *loopIndices = vrna_loopidx_from_ptable(pt);
move m3 ={3,7};
vrna_update_loop_indices(strlen(structure), loopIndices, &m3, pt);
int notEqual = 0;
for(int i =0; i < pt[0]; i++){
  if(loopIndices[i] != expectedLoopIndices3[i]){
    notEqual = 1;
  }
}
ck_assert_int_eq(notEqual, 0);

free (loopIndices);
free (pt);

//TODO: test all combinations of insertions and deletions

#test test_vrna_build_neighbors_iteratively

char *sequence =  "GGGAAACCCAACCUUU";
char *structure = ".(.....)........";
//                ".(.....).(.....)" after move
move currentMove = {10,16}; // for the second neighbor generation.
int expectedLength = 5;
int shifts = 0;
move expectedNeighbors[5] =
    { { 1, 9 }, { 3, 7 }, { -2, -8 }, { -10, -16 }, { 11, 15 } };
vrna_md_t md;
vrna_md_set_default (&md);
vrna_fold_compound_t *vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
short *previousStructure = vrna_ptable(structure);

move * previousMoves = vrna_build_neighbors(vc, previousStructure, shifts);
int length = 0;
for(move *m = previousMoves; m->bpLeft != 0; m++){
  length++;
}
int newLength;

move * neighbors = vrna_build_neighbors_iteratively(vc, &currentMove, previousStructure, previousMoves, length,
    &newLength, shifts);


//check if all neighbors are in the set.
int neighborsFound = 0;
for(int i = 0; i < expectedLength;i++){
  for(move *m = neighbors; m->bpLeft != 0; m++){
    if(expectedNeighbors[i].bpLeft == m->bpLeft && expectedNeighbors[i].bpRight == m->bpRight){
      neighborsFound++;
      break;
    }
  }
}

ck_assert_int_eq(neighborsFound, expectedLength);

vrna_fold_compound_free (vc);
free (previousStructure);
free(previousMoves);
free (neighbors);


//TODO: test more shift combinations

#test test_vrna_build_neighbors_iteratively_with_shifts

char *sequence =  "GAAACCAACCU";
char *structure = "(....).....";
//                "(....)(...)" after move
//                "123456789..
move currentMove = {7,11}; // for the second neighbor generation.
int expectedLength = 3;
int shifts = 1;
move expectedNeighbors[3] =
    { 
        { -1, -6 }, { -7, -11 }, 
        { 1, -5 }
    };
vrna_md_t md;
vrna_md_set_default (&md);
vrna_fold_compound_t *vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
short *previousStructure = vrna_ptable(structure);

int *loopIndices = vrna_loopidx_from_ptable(previousStructure);

move * previousMoves = vrna_build_neighbors(vc, previousStructure, shifts);
int length = 0;
for(move *m = previousMoves; m->bpLeft != 0; m++){
  length++;
}
int newLength;

move * neighbors = vrna_build_neighbors_iteratively(vc, &currentMove, previousStructure, previousMoves, length,
    &newLength, shifts);


//check if all neighbors are in the set.
int neighborsFound = 0;
for(int i = 0; i < expectedLength;i++){
  for(move *m = neighbors; m->bpLeft != 0; m++){
    if(expectedNeighbors[i].bpLeft == m->bpLeft && expectedNeighbors[i].bpRight == m->bpRight){
      neighborsFound++;
      break;
    }
  }
}

ck_assert_int_eq(neighborsFound, expectedLength);

vrna_fold_compound_free (vc);
free (previousStructure);
free (loopIndices);
free(previousMoves);
free (neighbors);
















