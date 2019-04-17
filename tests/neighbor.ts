#include <stdio.h>
#include <stdlib.h>
#include <ViennaRNA/landscape/neighbor.h>
#include <ViennaRNA/model.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/utils/structures.h>
#include <ViennaRNA/data_structures.h>
#include <stdarg.h>


/*************************
**  Utility Functions  **
*************************/

/**
 * @brief Test two arrays of C-style strings for deep equality. The results
 *  are reported via the testing framework (Check).
 * @param a First array of strings
 * @param size_a Size of first array
 * @param b Second array of strings
 * @param size_b Size of second array
 */
static void
deep_eq_str_ary(char  **a,
                int   size_a,
                char  **b,
                int   size_b)
{
  int i;

  ck_assert_int_eq(size_a, size_b);
  for (i = 0; i < size_a; i++)
    ck_assert_str_eq(a[i], b[i]);

  return;
}


/**
 * @brief Create dot-bracket structures from moves
 * @param moves list of vrna_move_t
 * @param length the number of moves.
 * @param pt_structure the RNA structure on which the moves can be applied.
 * @param seqLen the length of the RNA sequence.
 * @return Pointer to an array of structure dot-bracket strings
 */
static char **
createStructuresFromMoves(vrna_move_t *moves,
                          size_t      length,
                          const short *pt_structure,
                          int         seqLen)
{
  char  *cur      = NULL;
  char  **structs = (char **)malloc(3 * seqLen * seqLen * sizeof(char *));
  int   size      = 0;

  //while( cur = pop())
  for (vrna_move_t *m = moves; m->pos_5 != 0; m++) {
    short *pt_tmp = vrna_ptable_copy(pt_structure);
    vrna_move_apply(pt_tmp, m);
    cur = vrna_db_from_ptable(pt_tmp);

    structs[size++] = cur;
    free(pt_tmp);
  }

  return structs;
}


/************************************************************
* QuickSort stolen from
* http://www.comp.dit.ie/rlawlor/Alg_DS/sorting/quickSort.c
* Adapted to sort string arrays, core dumps fixed.
************************************************************/
static int
partition(char  *a[],
          int   l,
          int   r)
{
  int   i, j;
  char  *pivot, *t;

  pivot = a[l];
  i     = l;
  j     = r + 1;

  while (1) {
    /* FK: first compare i and r to prevent core dumps */
    do
      ++i;
    while (i <= r && strcmp(a[i], pivot) <= 0);
    do
      --j;
    while (strcmp(a[j], pivot) > 0);
    if (i >= j)
      break;

    t     = a[i];
    a[i]  = a[j];
    a[j]  = t;
  }
  t     = a[l];
  a[l]  = a[j];
  a[j]  = t;
  return j;
}


static void
quickSort(char  *a[],
          int   l,
          int   r)
{
  int j;

  if (l < r) {
    // divide and conquer
    j = partition(a, l, r);
    quickSort(a, l, j - 1);
    quickSort(a, j + 1, r);
  }
}


/**
 * @brief QuickSort an array of zero-terminated strings
 * @param ary Array of strings (i.e. array of char*)
 * @param size Number of strings in array
 */
static void
sort_str(char *ary[],
         int  size)
{
  quickSort(ary, 0, size - 1);
}


/**
 * @brief Test the function RNA2_move_it by comparing the generated neighbor
 *  to a passed list of structures. The passed structures are sorted, so their
 *  order is irrelevant. Results will be reported by the test framework (Check)
 * @param vc vrna_fold_compound_t
 * @param str Structure of which to generate neighbors
 * @param nb_size Number of passed neighbor structures (cf. '...')
 * @param ... several structure dot-bracket strings to which the generated
 *  neighbors should be compared to. This list will be sorted, so the order
 *  in which structures are passed is irrelevant.
 */
static void
test_RNA2_move_it(vrna_fold_compound_t  *vc,
                  char                  *str,
                  unsigned int          options,
                  int                   nb_size,
                  ...)
{
  int     i;
  va_list valist;   /* Variable args list */

  va_start(valist, nb_size);
  int     seqLen = vc->length;

  char    **gen_nb = NULL;  /* To store generated neighbors */
  int     gen_nb_size;      /* Number of generated structures */
  char    **nb = (char **)vrna_alloc(3 * seqLen * seqLen * sizeof(char *));

  for (i = 0; i < nb_size; i++)
    nb[i] = va_arg(valist, char *);
  sort_str(nb, nb_size);

  short       *pt_structure = vrna_ptable(str);
  vrna_move_t *neighbors    = vrna_neighbors(vc, pt_structure, options);
  int         length        = 0;
  for (vrna_move_t *m = neighbors; m->pos_5 != 0; m++)
    length++;

  gen_nb      = createStructuresFromMoves(neighbors, length, pt_structure, seqLen); /* Collect neighbors from stapel */
  gen_nb_size = length;
  sort_str(gen_nb, gen_nb_size);                                                    /* Sort generated neighbors */
  /*fprintf( stderr, "\"%s\" (input)\n", str); */
  /* for(i = 0; i<gen_nb_size; i++)      /* Output generated neighbors */
  /* fprintf( stderr, "\"%s\",\n", gen_nb[i]); */
  deep_eq_str_ary(gen_nb, gen_nb_size, nb, nb_size);

  for (int i = 0; i < length; i++)
    free(gen_nb[i]);
  free(gen_nb);
  free(nb);
  va_end(valist);
  free(pt_structure);
  vrna_move_list_free(neighbors);
}


#suite Neighbor

#test test_vrna_neighbors
{
  char                  *sequence             = "GGGAAACCCAACCUUU";
  char                  *structure            = ".(.....)........";
  int                   expectedLength        = 8;
  vrna_move_t           expectedNeighbors[8]  =
  { { 1, 9 }, { 3, 7 }, { -2, -8 }, { 10, 16 }, { 10, 15 }, { 10, 14 }, { 11, 16 }, { 11, 15 } };
  vrna_md_t             md;
  vrna_md_set_default(&md);
  vrna_fold_compound_t  *vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
  short                 *pt = vrna_ptable(structure);

  vrna_move_t           *neighbors = vrna_neighbors(vc, pt, VRNA_MOVESET_DEFAULT);

  //check if all neighbors are in the set.
  int                   neighborsFound = 0;
  for (int i = 0; i < expectedLength; i++) {
    for (vrna_move_t *m = neighbors; m->pos_5 != 0; m++) {
      if (expectedNeighbors[i].pos_5 == m->pos_5 && expectedNeighbors[i].pos_3 == m->pos_3) {
        neighborsFound++;
        break;
      }
    }
  }

  ck_assert_int_eq(neighborsFound, expectedLength);

  vrna_fold_compound_free(vc);
  free(pt);
  free(neighbors);
}

#test test_vrna_neighbors_with_shifts
{
  char                  *sequence   = "GAAACCAACCU";
  char                  *structure  = "(....).....";
  //                "(....)(...)" after vrna_move_t
  //                "123456789..
  int                   expectedLength        = 6;
  vrna_move_t           expectedNeighbors[6]  =
  {
    { -1, -6 }, { 7, 11 },
    { 1,  -5 }, { 1, -9 },{ 1, -10 }, { 1, -11 }     /* shifts */
  };
  vrna_md_t             md;
  vrna_md_set_default(&md);
  vrna_fold_compound_t  *vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
  short                 *pt = vrna_ptable(structure);

  vrna_move_t           *neighbors =
    vrna_neighbors(vc, pt, VRNA_MOVESET_DEFAULT | VRNA_MOVESET_SHIFT);

  //check if all neighbors are in the set.
  int                   neighborsFound = 0;
  for (int i = 0; i < expectedLength; i++) {
    for (vrna_move_t *m = neighbors; m->pos_5 != 0; m++) {
      if (expectedNeighbors[i].pos_5 == m->pos_5 && expectedNeighbors[i].pos_3 == m->pos_3) {
        neighborsFound++;
        break;
      }
    }
  }

  ck_assert_int_eq(neighborsFound, expectedLength);

  vrna_fold_compound_free(vc);
  free(pt);
  free(neighbors);
}


#test test_vrna_perform_move
{
  char        *sequence   = "GGGAAACCCAACCUUU";
  char        *structure  = "................";
  short       *pt         = vrna_ptable(structure);
  vrna_move_t m           = {
    2, 8
  };
  vrna_move_apply(pt, &m);
  char        *resultInsertion  = vrna_db_from_ptable(pt);
  char        *expectedResult   = ".(.....)........";
  //check if the vrna_move_t was inserted.
  ck_assert_str_eq(resultInsertion, expectedResult);

  //check if deletion is possible
  vrna_move_t m2 = {
    -2, -8
  };
  vrna_move_apply(pt, &m2);
  char        *resultDeletion = vrna_db_from_ptable(pt);
  ck_assert_str_eq(resultDeletion, structure);

  free(resultInsertion);
  free(resultDeletion);
  free(pt);
}


#test test_update_loop_indices_deletion
{
  char        *sequence   = "GGGAAACCCAACCUUU";
  char        *structure  = ".(.....)........";
  //test deletion
  vrna_move_t m = {
    -2, -8
  };
  int         expectedLoopIndices[17] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  };
  short       *pt           = vrna_ptable(structure);
  int         *loopIndices  = vrna_loopidx_from_ptable(pt);

  vrna_loopidx_update(loopIndices, pt, strlen(structure), &m);

  int         notEqual = 0;
  for (int i = 0; i < pt[0]; i++)
    if (loopIndices[i] != expectedLoopIndices[i])
      notEqual = 1;

  ck_assert_int_eq(notEqual, 0);

  free(loopIndices);
  free(pt);
}

#test test_update_loop_indices_insertion
{
  //test insertion
  //                "GGGAAACCCAACCUUU";
  char        *structure                = "................";
  int         expectedLoopIndices1[17]  = {
    1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0
  };
  short       *pt           = vrna_ptable(structure);
  int         *loopIndices  = vrna_loopidx_from_ptable(pt);
  vrna_move_t m1            = {
    2, 8
  };
  vrna_loopidx_update(loopIndices, pt, strlen(structure), &m1);
  int         notEqual = 0;
  for (int i = 0; i < pt[0]; i++)
    if (loopIndices[i] != expectedLoopIndices1[i])
      notEqual = 1;

  ck_assert_int_eq(notEqual, 0);

  free(loopIndices);
  free(pt);
}


#test test_update_loop_indices_insertion_after_insertion
{
  //test insertion after insertion
  //          "GGGAAACCCAACCUUU";
  char        *structure = ".(.....)........";
  //          ".(.....)..(...)." after vrna_move_t
  int         expectedLoopIndices2[17] = {
    2, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 2, 2, 2, 2, 2, 0
  };
  short       *pt           = vrna_ptable(structure);
  int         *loopIndices  = vrna_loopidx_from_ptable(pt);
  vrna_move_t m2            = {
    11, 15
  };
  vrna_loopidx_update(loopIndices, pt, strlen(structure), &m2);
  int         notEqual = 0;
  for (int i = 0; i < pt[0]; i++)
    if (loopIndices[i] != expectedLoopIndices2[i])
      notEqual = 1;

  ck_assert_int_eq(notEqual, 0);

  free(loopIndices);
  free(pt);
}


#test test_update_loop_indices_insertion_within_insertion
{
  //test insertion within insertion
  //          "GGGAAACCCAACCUUU";
  char        *structure = ".(.....)........";
  //          ".((...))........" after vrna_move_t
  int         expectedLoopIndices3[17] = {
    2, 0, 1, 2, 2, 2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0
  };
  short       *pt           = vrna_ptable(structure);
  int         *loopIndices  = vrna_loopidx_from_ptable(pt);
  vrna_move_t m3            = {
    3, 7
  };
  vrna_loopidx_update(loopIndices, pt, strlen(structure), &m3);
  int         notEqual = 0;
  for (int i = 0; i < pt[0]; i++)
    if (loopIndices[i] != expectedLoopIndices3[i])
      notEqual = 1;

  ck_assert_int_eq(notEqual, 0);

  free(loopIndices);
  free(pt);
}

//TODO: test all combinations of insertions and deletions

#test test_vrna_neighbors_successive
{
  char                  *sequence   = "GGGAAACCCAACCUUU";
  char                  *structure  = ".(.....)........";
  //                ".(.....).(.....)" after vrna_move_t
  vrna_move_t           currentMove = {
    10, 16
  };                               // for the second neighbor generation.
  int                   expectedLength        = 5;
  vrna_move_t           expectedNeighbors[5]  =
  { { 1, 9 }, { 3, 7 }, { -2, -8 }, { -10, -16 }, { 11, 15 } };
  vrna_md_t             md;
  vrna_md_set_default(&md);
  vrna_fold_compound_t  *vc =
    vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
  short                 *previousStructure = vrna_ptable(structure);

  vrna_move_t           *previousMoves =
    vrna_neighbors(vc, previousStructure, VRNA_MOVESET_DEFAULT);
  int                   length = 0;
  for (vrna_move_t *m = previousMoves; m->pos_5 != 0; m++)
    length++;
  int                   newLength;

  vrna_move_t           *neighbors = vrna_neighbors_successive(vc,
                                                               &currentMove,
                                                               previousStructure,
                                                               previousMoves,
                                                               length,
                                                               &newLength,
                                                               VRNA_MOVESET_DEFAULT);


  //check if all neighbors are in the set.
  int neighborsFound = 0;
  for (int i = 0; i < expectedLength; i++) {
    for (vrna_move_t *m = neighbors; m->pos_5 != 0; m++) {
      if (expectedNeighbors[i].pos_5 == m->pos_5 && expectedNeighbors[i].pos_3 == m->pos_3) {
        neighborsFound++;
        break;
      }
    }
  }

  ck_assert_int_eq(neighborsFound, expectedLength);

  vrna_fold_compound_free(vc);
  free(previousStructure);
  free(previousMoves);
  free(neighbors);
}

#test test_vrna_neighbors_successive_deletions_only
{
  char                  *sequence   = "GGGAAACCCAACCUUU";
  char                  *structure  = ".(.....)........";
  //                ".(.....).(.....)" after vrna_move_t
  vrna_move_t           currentMove = {
    10, 16
  };                               // for the second neighbor generation.
  int                   expectedLength        = 2;
  vrna_move_t           expectedNeighbors[2]  =
  { { -2, -8 }, { -10, -16 } };
  vrna_md_t             md;
  vrna_md_set_default(&md);
  vrna_fold_compound_t  *vc =
    vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
  short                 *previousStructure = vrna_ptable(structure);

  vrna_move_t           *previousMoves =
    vrna_neighbors(vc, previousStructure, VRNA_MOVESET_DELETION);
  int                   length = 0;
  for (vrna_move_t *m = previousMoves; m->pos_5 != 0; m++)
    length++;
  int                   newLength;

  vrna_move_t           *neighbors = vrna_neighbors_successive(vc,
                                                               &currentMove,
                                                               previousStructure,
                                                               previousMoves,
                                                               length,
                                                               &newLength,
                                                               VRNA_MOVESET_DELETION);


  //check if all neighbors are in the set.
  int neighborsFound = 0;
  for (int i = 0; i < expectedLength; i++) {
    for (vrna_move_t *m = neighbors; m->pos_5 != 0; m++) {
      if (expectedNeighbors[i].pos_5 == m->pos_5 && expectedNeighbors[i].pos_3 == m->pos_3) {
        neighborsFound++;
        break;
      }
    }
  }

  ck_assert_int_eq(neighborsFound, expectedLength);

  vrna_fold_compound_free(vc);
  free(previousStructure);
  free(previousMoves);
  free(neighbors);
}

#test test_vrna_neighbors_successive_insertions_only
{
  char                  *sequence   = "GGGAAACCCAACCUUU";
  char                  *structure  = ".(.....)........";
  //                ".(.....).(.....)" after vrna_move_t
  vrna_move_t           currentMove = {
    10, 16
  };                               // for the second neighbor generation.
  int                   expectedLength        = 3;
  vrna_move_t           expectedNeighbors[3]  =
  { { 1, 9 }, { 3, 7 }, { 11, 15 } };
  vrna_md_t             md;
  vrna_md_set_default(&md);
  vrna_fold_compound_t  *vc =
    vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
  short                 *previousStructure = vrna_ptable(structure);

  vrna_move_t           *previousMoves =
    vrna_neighbors(vc, previousStructure, VRNA_MOVESET_DEFAULT);
  int                   length = 0;
  for (vrna_move_t *m = previousMoves; m->pos_5 != 0; m++)
    length++;
  int                   newLength;

  vrna_move_t           *neighbors = vrna_neighbors_successive(vc,
                                                               &currentMove,
                                                               previousStructure,
                                                               previousMoves,
                                                               length,
                                                               &newLength,
                                                               VRNA_MOVESET_DEFAULT);


  //check if all neighbors are in the set.
  int neighborsFound = 0;
  for (int i = 0; i < expectedLength; i++) {
    for (vrna_move_t *m = neighbors; m->pos_5 != 0; m++) {
      if (expectedNeighbors[i].pos_5 == m->pos_5 && expectedNeighbors[i].pos_3 == m->pos_3) {
        neighborsFound++;
        break;
      }
    }
  }

  ck_assert_int_eq(neighborsFound, expectedLength);

  vrna_fold_compound_free(vc);
  free(previousStructure);
  free(previousMoves);
  free(neighbors);
}

//TODO: test more shift combinations

#test test_vrna_neighbors_successive_with_shifts
{
  char                  *sequence   = "GAAACCAACCU";
  char                  *structure  = "(....).....";
  //                "(....)(...)" after vrna_move_t
  //                "123456789..
  vrna_move_t           currentMove = {
    7, 11
  };                              // for the second neighbor generation.
  int                   expectedLength        = 3;
  vrna_move_t           expectedNeighbors[3]  =
  {
    { -1, -6 }, { -7, -11 },
    { 1,  -5 }
  };
  vrna_md_t             md;
  vrna_md_set_default(&md);
  vrna_fold_compound_t  *vc =
    vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);
  short                 *previousStructure = vrna_ptable(structure);

  vrna_move_t           *previousMoves = vrna_neighbors(vc,
                                                        previousStructure,
                                                        VRNA_MOVESET_DEFAULT | VRNA_MOVESET_SHIFT);
  int                   length = 0;
  for (vrna_move_t *m = previousMoves; m->pos_5 != 0; m++)
    length++;
  int                   newLength;

  vrna_move_t           *neighbors = vrna_neighbors_successive(vc,
                                                               &currentMove,
                                                               previousStructure,
                                                               previousMoves,
                                                               length,
                                                               &newLength,
                                                               VRNA_MOVESET_DEFAULT |
                                                               VRNA_MOVESET_SHIFT);


  //check if all neighbors are in the set.
  int neighborsFound = 0;
  for (int i = 0; i < expectedLength; i++) {
    for (vrna_move_t *m = neighbors; m->pos_5 != 0; m++) {
      if (expectedNeighbors[i].pos_5 == m->pos_5 && expectedNeighbors[i].pos_3 == m->pos_3) {
        neighborsFound++;
        break;
      }
    }
  }

  ck_assert_int_eq(neighborsFound, expectedLength);

  vrna_fold_compound_free(vc);
  free(previousStructure);
  free(previousMoves);
  free(neighbors);
}

/* Test generation of neighbor structures, without shift moves, with lonely pairs */
#test test_rnamoves_noshift_lp
{
  //int noLP = 0;
  //int shift = 0;
  unsigned int          options = VRNA_MOVESET_DEFAULT;

  char                  *seq    = "GGGGUUUCCCC";
  char                  *str    = "((((...))))";
  int                   nb_size = 4;


  vrna_md_t             md;
  vrna_md_set_default(&md);
  vrna_fold_compound_t  *vc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);

  test_RNA2_move_it(vc, str, options, nb_size,
                    "(((.....)))",
                    "((.(...).))",
                    "(.((...)).)",
                    ".(((...)))."
                    );

  vrna_fold_compound_free(vc);

  seq     = "GGUUUCC";
  str     = ".......";
  nb_size = 5;                  /* Number of passed neighbors */

  vc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);
  test_RNA2_move_it(vc, str, options, nb_size,
                    ".(....)",
                    ".(...).",
                    "(.....)",
                    "(....).",
                    "(...).."
                    );

  vrna_fold_compound_free(vc);
}


/* Test generation of neighbor structures, with shift moves, with lonely pairs */
#test test_rnamoves_shift_lp
{
  //int shift = 1;
  //int noLP = 0;
  unsigned int          options = VRNA_MOVESET_DEFAULT | VRNA_MOVESET_SHIFT;

  char                  *seq  = "GGGGGGGGGGGGAAAUUUUUUUGUUUU";
  char                  *str  = ".(((....(((.....)))....))).";

  vrna_md_t             md;
  vrna_md_set_default(&md);
  vrna_fold_compound_t  *vc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);

  int                   nb_size = 40;
  test_RNA2_move_it(vc, str, options, nb_size,
                    ".(((...((((.....))))...))).", ".(((.(..(((.....)))..).))).",
                    "((((....(((.....)))....))))", ".(((....((((...))))....))).",
                    "..((....(((.....)))....))..", ".(((.....((.....)).....))).",
                    ".(((....((.......))....))).", ".(((..(.(((.....))).)..))).",
                    ".((.....(((.....))).....)).", ".(.(....(((.....)))....).).",
                    ".(((....(.(.....).)....))).", ".(((....(((....).))....))).",
                    ".(((.(..(((.....))).)..))).", ".(((..(.(((.....)))..).))).",
                    ".(((....((.(....)))....))).", ".(((..(.(((.....))))...))).",
                    ".(((...((((.....))).)..))).", ".(((....(((.....)))....)).)",
                    ".((.(...(((.....)))....))).", ".(((...(.((.....)))....))).",
                    ".(((....(((.....)).)...))).", "(.((....(((.....)))....))).",
                    ".((((...(((.....)))..).))).", ".((((...(((.....))).)..))).",
                    ".(((.(..(((.....))))...))).", ".(((...((((.....)))..).))).",
                    ".(((....((..(...)))....))).", ".((..(..(((.....)))....))).",
                    ".(((..(..((.....)))....))).", ".(((....(((.....))..)..))).",
                    ".((((...(((.....))))...))).", ".(((....(((.....)))..)..)).",
                    ".((....((((.....)))....))).", ".(((....(((.....))))....)).",
                    ".((((....((.....)))....))).", ".((...(.(((.....)))....))).",
                    ".(((....(((.....))).)...)).", ".(((.(...((.....)))....))).",
                    ".(((....(((.....))...).))).", ".(((.....((.....))(...))))."
                    );
  vrna_fold_compound_free(vc);
}

/* Test generation of neighbor structures, without shift moves, without lonely pairs */
#test test_rnamoves_noshift_nolp
{
  //int shift = 0;
  //int noLP = 1;
  unsigned int          options = VRNA_MOVESET_DEFAULT | VRNA_MOVESET_NO_LP;

  char                  *seq  = "GGGGUUUCCCC";
  char                  *str  = "((((...))))";

  vrna_md_t             md;
  vrna_md_set_default(&md);
  vrna_fold_compound_t  *vc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);

  int                   nb_size = 2;
  test_RNA2_move_it(vc, str, options, nb_size,
                    "(((.....)))",
                    ".(((...)))."
                    );

  str     = "...........";
  nb_size = 12;
  test_RNA2_move_it(vc, str, options, nb_size,
                    "((......)).", ".((....))..", ".((.....)).",
                    "((.....))..", "..((....)).", ".((......))",
                    "..((...))..", "((.......))", "..((.....))",
                    "((....))...", ".((...))...", "((...))...."
                    );
  vrna_fold_compound_free(vc);
}

/* Test generation of neighbor structures, with shift moves, without lonely pairs */
#test test_rnamoves_shift_nolp
{
  //int shift = 1;
  //int noLP = 1;
  unsigned int          options = VRNA_MOVESET_DEFAULT | VRNA_MOVESET_SHIFT | VRNA_MOVESET_NO_LP;

  char                  *seq  = "GGGGGGGGUUUUUUU";
  char                  *str  = "((..(((...)))))";
  //               01234567890123456789012
  //               0         1         2

  vrna_md_t             md;
  vrna_md_set_default(&md);
  vrna_fold_compound_t  *vc = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);

  int                   nb_size = 4;
  test_RNA2_move_it(vc, str, options, nb_size,
                    "....(((...)))..",
                    "((..((.....))))",
                    "(((..((...)))))",
                    "((...((...)).))"
                    );

  vrna_fold_compound_free(vc);

  seq     = "AAGUGAUACCAGCAUCGUCUUGAUGCCCUUGGCAGCACUUCAU";
  str     = "(((((((......)))).)))((((((...)))).....))..";
  nb_size = 9;
  vc      = vrna_fold_compound(seq, &md, VRNA_OPTION_EVAL_ONLY);
  test_RNA2_move_it(vc, str, options, nb_size,
                    "(((((((......))))).))((((((...)))).....))..",
                    "(((((((......)))).)))(((((.....))).....))..",
                    "(((((((......)))).)))((.(((...)))......))..",
                    "(((((((......)))).)))..((((...)))).........",
                    "((((((........))).)))((((((...)))).....))..",
                    "(((.(((......)))..)))((((((...)))).....))..",
                    "((.((((......))))..))((((((...)))).....))..",
                    ".((((((......)))).))(((((((...)))).....))).",
                    ".((((((......)))).)).((((((...)))).....)).."
                    );
  vrna_fold_compound_free(vc);
}
