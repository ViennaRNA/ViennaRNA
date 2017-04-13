#include <stdio.h>
#include <stdlib.h>
#include <ViennaRNA/neighbor.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/structure_utils.h>
#include <ViennaRNA/data_structures.h>

#suite Neighbor

#test test_vrna_neighbors

char *sequence = "GGGAAACCCAACCUUU";
char *structure = ".(.....)........";
int expectedLength = 8;
vrna_move_t expectedNeighbors[8] =
    { { 1, 9 }, { 3, 7 }, { -2, -8 }, { 10, 16 }, { 10, 15 }, { 10, 14 }, { 11, 16 }, { 11, 15 } };
vrna_md_t md;
vrna_md_set_default (&md);
vrna_fold_compound_t *vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
short *pt = vrna_ptable(structure);

vrna_move_t * neighbors = vrna_neighbors(vc, pt, VRNA_MOVESET_DEFAULT);

//check if all neighbors are in the set.
int neighborsFound = 0;
for(int i = 0; i < expectedLength;i++){
  for(vrna_move_t *m = neighbors; m->pos_5 != 0; m++){
    if(expectedNeighbors[i].pos_5 == m->pos_5 && expectedNeighbors[i].pos_3 == m->pos_3){
      neighborsFound++;
      break;
    }
  }
}

ck_assert_int_eq(neighborsFound, expectedLength);

vrna_fold_compound_free (vc);
free (pt);
free (neighbors);


#test test_vrna_neighbors_with_shifts

char *sequence =  "GAAACCAACCU";
char *structure = "(....).....";
//                "(....)(...)" after vrna_move_t
//                "123456789..
int expectedLength = 6;
vrna_move_t expectedNeighbors[6] =
    { 
        { -1, -6 }, { 7, 11 }, 
        { 1, -5 }, { 1, -9 }, { 1, -10 }, { 1, -11 } /* shifts */
    };
vrna_md_t md;
vrna_md_set_default (&md);
vrna_fold_compound_t *vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
short *pt = vrna_ptable(structure);

vrna_move_t * neighbors = vrna_neighbors(vc, pt, VRNA_MOVESET_DEFAULT | VRNA_MOVESET_SHIFT);

//check if all neighbors are in the set.
int neighborsFound = 0;
for(int i = 0; i < expectedLength;i++){
  for(vrna_move_t *m = neighbors; m->pos_5 != 0; m++){
    if(expectedNeighbors[i].pos_5 == m->pos_5 && expectedNeighbors[i].pos_3 == m->pos_3){
      neighborsFound++;
      break;
    }
  }
}

ck_assert_int_eq(neighborsFound, expectedLength);

vrna_fold_compound_free (vc);
free (pt);
free (neighbors);


#test test_vrna_perform_move

char *sequence = "GGGAAACCCAACCUUU";
char *structure = "................";
short *pt = vrna_ptable(structure);
vrna_move_t m = {2, 8};
vrna_move_apply(pt, &m);
char * resultInsertion = vrna_db_from_ptable(pt);
char *expectedResult = ".(.....)........";
//check if the vrna_move_t was inserted.
ck_assert_str_eq( resultInsertion, expectedResult );

//check if deletion is possible
vrna_move_t m2 = {-2, -8};
vrna_move_apply(pt, &m2);
char * resultDeletion = vrna_db_from_ptable(pt);
ck_assert_str_eq( resultDeletion, structure );

free (resultInsertion);
free (resultDeletion);


#test test_update_loop_indices_deletion

char *sequence = "GGGAAACCCAACCUUU";
char *structure = ".(.....)........";
//test deletion
vrna_move_t m = { -2, -8 };
int expectedLoopIndices[17] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
short *pt = vrna_ptable(structure);
int *loopIndices = vrna_loopidx_from_ptable(pt);

vrna_loopidx_update(loopIndices, pt, strlen(structure), &m);

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
vrna_move_t m1 ={2,8};
vrna_loopidx_update(loopIndices, pt, strlen(structure), &m1);
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
//          ".(.....)..(...)." after vrna_move_t
int expectedLoopIndices2[17] = { 2, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 2, 2, 2, 2, 2, 0 };
short *pt = vrna_ptable(structure);
int *loopIndices = vrna_loopidx_from_ptable(pt);
vrna_move_t m2 ={11,15};
vrna_loopidx_update(loopIndices, pt, strlen(structure), &m2);
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
//          ".((...))........" after vrna_move_t
int expectedLoopIndices3[17] = { 2, 0, 1, 2, 2, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0 };
short *pt = vrna_ptable(structure);
int *loopIndices = vrna_loopidx_from_ptable(pt);
vrna_move_t m3 ={3,7};
vrna_loopidx_update(loopIndices, pt, strlen(structure), &m3);
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

#test test_vrna_neighbors_successive

char *sequence =  "GGGAAACCCAACCUUU";
char *structure = ".(.....)........";
//                ".(.....).(.....)" after vrna_move_t
vrna_move_t currentMove = {10,16}; // for the second neighbor generation.
int expectedLength = 5;
vrna_move_t expectedNeighbors[5] =
    { { 1, 9 }, { 3, 7 }, { -2, -8 }, { -10, -16 }, { 11, 15 } };
vrna_md_t md;
vrna_md_set_default (&md);
vrna_fold_compound_t *vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
short *previousStructure = vrna_ptable(structure);

vrna_move_t * previousMoves = vrna_neighbors(vc, previousStructure, VRNA_MOVESET_DEFAULT);
int length = 0;
for(vrna_move_t *m = previousMoves; m->pos_5 != 0; m++){
  length++;
}
int newLength;

vrna_move_t * neighbors = vrna_neighbors_successive(vc, &currentMove, previousStructure, previousMoves, length,
    &newLength, VRNA_MOVESET_DEFAULT);


//check if all neighbors are in the set.
int neighborsFound = 0;
for(int i = 0; i < expectedLength;i++){
  for(vrna_move_t *m = neighbors; m->pos_5 != 0; m++){
    if(expectedNeighbors[i].pos_5 == m->pos_5 && expectedNeighbors[i].pos_3 == m->pos_3){
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


#test test_vrna_neighbors_successive_deletions_only

char *sequence =  "GGGAAACCCAACCUUU";
char *structure = ".(.....)........";
//                ".(.....).(.....)" after vrna_move_t
vrna_move_t currentMove = {10,16}; // for the second neighbor generation.
int expectedLength = 2;
vrna_move_t expectedNeighbors[2] =
    { { -2, -8 }, { -10, -16 }};
vrna_md_t md;
vrna_md_set_default (&md);
vrna_fold_compound_t *vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
short *previousStructure = vrna_ptable(structure);

vrna_move_t * previousMoves = vrna_neighbors(vc, previousStructure, VRNA_MOVESET_DELETION);
int length = 0;
for(vrna_move_t *m = previousMoves; m->pos_5 != 0; m++){
  length++;
}
int newLength;

vrna_move_t * neighbors = vrna_neighbors_successive(vc, &currentMove, previousStructure, previousMoves, length,
    &newLength, VRNA_MOVESET_DELETION);


//check if all neighbors are in the set.
int neighborsFound = 0;
for(int i = 0; i < expectedLength;i++){
  for(vrna_move_t *m = neighbors; m->pos_5 != 0; m++){
    if(expectedNeighbors[i].pos_5 == m->pos_5 && expectedNeighbors[i].pos_3 == m->pos_3){
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


#test test_vrna_neighbors_successive_insertions_only

char *sequence =  "GGGAAACCCAACCUUU";
char *structure = ".(.....)........";
//                ".(.....).(.....)" after vrna_move_t
vrna_move_t currentMove = {10,16}; // for the second neighbor generation.
int expectedLength = 3;
vrna_move_t expectedNeighbors[3] =
    { { 1, 9 }, { 3, 7 }, { 11, 15 } };
vrna_md_t md;
vrna_md_set_default (&md);
vrna_fold_compound_t *vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
short *previousStructure = vrna_ptable(structure);

vrna_move_t * previousMoves = vrna_neighbors(vc, previousStructure, VRNA_MOVESET_DEFAULT);
int length = 0;
for(vrna_move_t *m = previousMoves; m->pos_5 != 0; m++){
  length++;
}
int newLength;

vrna_move_t * neighbors = vrna_neighbors_successive(vc, &currentMove, previousStructure, previousMoves, length,
    &newLength, VRNA_MOVESET_DEFAULT);


//check if all neighbors are in the set.
int neighborsFound = 0;
for(int i = 0; i < expectedLength;i++){
  for(vrna_move_t *m = neighbors; m->pos_5 != 0; m++){
    if(expectedNeighbors[i].pos_5 == m->pos_5 && expectedNeighbors[i].pos_3 == m->pos_3){
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

#test test_vrna_neighbors_successive_with_shifts

char *sequence =  "GAAACCAACCU";
char *structure = "(....).....";
//                "(....)(...)" after vrna_move_t
//                "123456789..
vrna_move_t currentMove = {7,11}; // for the second neighbor generation.
int expectedLength = 3;
vrna_move_t expectedNeighbors[3] =
    { 
        { -1, -6 }, { -7, -11 }, 
        { 1, -5 }
    };
vrna_md_t md;
vrna_md_set_default (&md);
vrna_fold_compound_t *vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
short *previousStructure = vrna_ptable(structure);

vrna_move_t * previousMoves = vrna_neighbors(vc, previousStructure, VRNA_MOVESET_DEFAULT | VRNA_MOVESET_SHIFT);
int length = 0;
for(vrna_move_t *m = previousMoves; m->pos_5 != 0; m++){
  length++;
}
int newLength;

vrna_move_t * neighbors = vrna_neighbors_successive(vc, &currentMove, previousStructure, previousMoves, length,
    &newLength, VRNA_MOVESET_DEFAULT | VRNA_MOVESET_SHIFT);


//check if all neighbors are in the set.
int neighborsFound = 0;
for(int i = 0; i < expectedLength;i++){
  for(vrna_move_t *m = neighbors; m->pos_5 != 0; m++){
    if(expectedNeighbors[i].pos_5 == m->pos_5 && expectedNeighbors[i].pos_3 == m->pos_3){
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
















