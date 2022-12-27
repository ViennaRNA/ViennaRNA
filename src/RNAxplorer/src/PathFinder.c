/*
 * PathFinder algorithm
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/2Dfold.h>
#include <ViennaRNA/walk.h>

#include "RNAwalk.h"
#include "PathFinder.h"


typedef char *(get_transient_structure_func)(vrna_fold_compound_t    *fc,
                                              const char              *start_structure,
                                              rnax_path_finder_opt_t  *options);

static unsigned int
find_saddle_point(vrna_path_t *folding_path);


static vrna_path_t *
concat_path(vrna_path_t *target,
            vrna_path_t *extension);


static vrna_path_t *
clone_path_element(vrna_path_t *source);


static char *
get_transient_gradient_walk(vrna_fold_compound_t    *fc,
                            const char              *start_structure,
                            rnax_path_finder_opt_t  *options);


static char *
get_transient_monte_carlo(vrna_fold_compound_t    *fc,
                          const char              *start_structure,
                          rnax_path_finder_opt_t  *options);


rnax_path_finder_opt_t *
rnax_path_finder_options(void)
{
  rnax_path_finder_opt_t *options = (rnax_path_finder_opt_t *)vrna_alloc(sizeof(rnax_path_finder_opt_t));

  options->method       = GRADIENT_WALK;
  options->iterations   = 1;
  options->max_keep     = 10;
  options->storage_size = 10;
  options->max_d1       = 5;
  options->max_d2       = 5;
  options->max_paths    = 1;

  vrna_md_set_default(&(options->md));

  return options;
}


vrna_path_t *
rnax_path_finder( const char              *seq,
                  const char              *s_source,
                  const char              *s_target,
                  rnax_path_finder_opt_t  *opt)
{
  char                          *transient_structure;
  unsigned int                  n, i, saddle_pos, s1_pos, s2_pos;
  double                        saddle_en, s1_en, s2_en, transient_en;
  vrna_path_t                   *refolding_path, *p1, *p2, **alternative_paths;
  rnax_path_finder_opt_t        *options;
  get_transient_structure_func  *get_transient_structures;

  refolding_path = NULL;

  if ((seq) && (s_source) && (s_target)) {
    transient_structure = NULL;
    p1                   = NULL;
    p2                   = NULL;
    n                    = strlen(seq);

    if ((strlen(s_source) != n) || (strlen(s_target) != n)) {
      vrna_message_warning("rnax_path_finder: lengths of structures do not match length of sequence!");
      return NULL;
    }

    options = (opt) ? opt : rnax_path_finder_options();

    switch (options->method) {
      case RNAX_PATHFINDER_SADDLE_GRADIENT_WALK:
        get_transient_structures = &get_transient_gradient_walk;
        break;

      case RNAX_PATHFINDER_SADDLE_MONTE_CARLO: /* fall through */
      case RNAX_PATHFINDER_SADDLE_MONTE_CARLO_SA:
        get_transient_structures = &get_transient_monte_carlo;
        break;

      default:
        break;
    }

    initRNAWalk (seq,  &(options->md));

    vrna_fold_compound_t *fc = vrna_fold_compound(seq, &(options->md), VRNA_OPTION_DEFAULT | VRNA_OPTION_EVAL_ONLY);

    /* 1st step, compute direct folding path */
    refolding_path  = vrna_path_findpath(fc, s_source, s_target, options->max_keep);
    saddle_pos      = find_saddle_point(refolding_path);
    saddle_en       = refolding_path[saddle_pos].en;

    /* 2nd step, collect alternative routes */
    alternative_paths     = (vrna_path_t **)vrna_alloc(sizeof(vrna_path_t *) * (options->storage_size));
    alternative_paths[0]  = vrna_path_findpath(fc, s_source, s_target, options->max_keep);


    for (i = 0; i < options->iterations; i++) {
      /* obtain new mesh point for detour path by starting a walk from current saddle */
      transient_structure = get_transient_structures(fc,
                                                      refolding_path[saddle_pos].s,
                                                      options);

      transient_en = vrna_eval_structure(fc, transient_structure);

      if ((transient_en < saddle_en) && (strcmp(refolding_path[saddle_pos].s, transient_structure) != 0)) {
        /* compute direct paths to new mesh point */
        p1 = vrna_path_findpath_ub(fc, s_source, transient_structure, options->max_keep, saddle_en);
        p2 = vrna_path_findpath_ub(fc, transient_structure, s_target, options->max_keep, saddle_en);

        if ((p1) && (p2)) {
          /* accept new path if saddle is lower than before */
          s1_pos = find_saddle_point(p1);
          s2_pos = find_saddle_point(p2);
          s1_en = p1[s1_pos].en;
          s2_en = p2[s2_pos].en;
          if ((s1_en < saddle_en) && (s2_en < saddle_en)) {
            /* release memory of previous path */
            free_path(refolding_path);

            refolding_path = concat_path(p1, p2 + 1);
            free(p2);

            p1          = NULL;
            p2          = NULL;
            saddle_pos  = find_saddle_point(refolding_path);
            saddle_en   = refolding_path[saddle_pos].en;
          } else {
            break;
          }
        } else {
          break;
        }
      } else {
        break;
      }

      free(transient_structure);
      free(p1);
      free(p2);

      transient_structure = NULL;
      p1                  = NULL;
      p2                  = NULL;
    }

    free(transient_structure);
    free(p1);
    free(p2);
    vrna_fold_compound_free(fc);

    if (!opt)
      free(options);
  }

  return refolding_path;
}


static char *
get_transient_gradient_walk(vrna_fold_compound_t    *fc,
                            const char              *start_structure,
                            rnax_path_finder_opt_t  *options)
{
  char *structure;
  short *pt       = vrna_ptable(start_structure);
  (void)vrna_path(fc, pt, 0, VRNA_PATH_DEFAULT | VRNA_PATH_NO_TRANSITION_OUTPUT);
  structure   = vrna_db_from_ptable(pt);
  free(pt);

  return structure;
}


static char *
get_transient_monte_carlo(vrna_fold_compound_t    *fc,
                          const char              *start_structure,
                          rnax_path_finder_opt_t  *options)
{
  char *structure;
    structure = structureWalk(fc->sequence,
                                  start_structure,
                                  options->method);

  return structure;
}


/**
 * Computes the direct folding path between the given structures. This path is used as template for an alternative path search.
 * The alternative path consists of meshpoints. The first meshpoint is the saddlepoint on the direct path. All other meshpoints are
 * computed depending on the number of maxIterations.
 * @param seq - the RNA sequence (ACGU)
 * @param s1 - first structure in dot-bracket format.
 * @param s2 - second structure in dot-bracket format.
 * @param maxIterations - maximal iterations of structurewalks.
 * @param maxKeep - how many structures are being kept (get_path in findpath.c)
 * @param method - gradient or random walk (GRADIENT_WALK, MC_METROPOLIS)
 * @param maxStorage - maximal number of meshpoints (see insert_meshpoint in meshpoint.c)
 */
void
levelSaddlePoint(const char *seq, const char *s1, const char *s2, int maxIterations, int maxKeep, int method, int maxStorage)
{

  int iterator = maxIterations;
  vrna_path_t *foldingPath = get_path (seq, s1, s2, maxKeep/*, &steps, circ*/);
  vrna_path_t *Saddle = getSaddlePoint (foldingPath/*, steps*/);

  float oldSaddleEn = Saddle->en;
  fprintf (stdout, "old Path:\nbarrier: %6.2f\n\n", Saddle->en);
  int d;
  for (d = 0; foldingPath[d].s; d++) {
    fprintf (stdout, "%s %6.2f\n", foldingPath[d].s, foldingPath[d].en);
  }
  vrna_path_t *newLeftSaddle = NULL;
  vrna_path_t *newRightSaddle = NULL;
  vrna_path_t *path_left, *path_right;
  float newSaddleEn = oldSaddleEn;
  meshpoint_list bestMeshPoints;
  init_meshpoint_list (&bestMeshPoints);
  fprintf (stdout, "\nsearching for alternative paths...");
  fflush (stdout);

  vrna_md_t md;
  vrna_md_set_default (&md);
  /* set user-defined model details */
  md.circ = circ;
  md.uniq_ML = 1;
  initRNAWalk (seq, &md);
  while (iterator > 0) {
    char *newSaddle = structureWalk (seq, Saddle->s, method);
    if (strcmp (s1, newSaddle) != 0) {
      path_left = get_path (seq, s1, newSaddle, maxKeep/*, &steps1, circ*/);
      path_right = get_path (seq, newSaddle, s2, maxKeep/*, &steps2, circ*/);

      if (newLeftSaddle != NULL) {
        if (newLeftSaddle->s)
          free (newLeftSaddle->s);
        free (newLeftSaddle);
      }
      if (newRightSaddle != NULL) {
        if (newRightSaddle->s)
          free (newRightSaddle->s);
        free (newRightSaddle);
      }
      newLeftSaddle = getSaddlePoint (path_left/*, steps1*/);
      newRightSaddle = getSaddlePoint (path_right/*, steps2*/);
      newSaddleEn = MAX2(newLeftSaddle->en, newRightSaddle->en);
    }
    else {
      continue;
    }
    iterator--;

    if (newSaddleEn < oldSaddleEn) {
      //TODO: -GE: is newSaddle the correct meshpoint? Or should it be newLeftSaddle or newRightSaddle?
      insert_meshpoint (newSaddle, newSaddleEn, &bestMeshPoints, maxStorage);
    }
    else {
      free (newSaddle);
    }
    free_path (path_left);
    free_path (path_right);
    fprintf (stdout, "\rsearching for alternative paths...(%2.1f %% done)", (float) (maxIterations - iterator) / (float) maxIterations * 100.);
    fflush (stdout);
  }

  if (bestMeshPoints.count > 0) {
    path_left = get_path (seq, s1, (bestMeshPoints.first)->s, maxKeep/*, &steps1, circ*/);
    path_right = get_path (seq, (bestMeshPoints.first)->s, s2, maxKeep/*, &steps2, circ*/);

    if (newLeftSaddle != NULL) {
      if (newLeftSaddle->s)
        free (newLeftSaddle->s);
      free (newLeftSaddle);
    }
    if (newRightSaddle != NULL) {
      if (newRightSaddle->s)
        free (newRightSaddle->s);
      free (newRightSaddle);
    }
    newLeftSaddle = getSaddlePoint (path_left/*, steps1*/);
    newRightSaddle = getSaddlePoint (path_right/*, steps2*/);
    newSaddleEn = MAX2(newLeftSaddle->en, newRightSaddle->en);

    fprintf (stdout, "\nnewPath with barrier: %6.2f\n\n", newSaddleEn);
    for (d = 0; path_left[d].s; d++) {
      fprintf (stdout, "%s %6.2f\n", path_left[d].s, path_left[d].en);
    }
    for (d = 1; path_right[d].s; d++) {
      fprintf (stdout, "%s %6.2f\n", path_right[d].s, path_right[d].en);
    }
  }
  else {
    fprintf (stdout, "no better path found...\n :-/\n");
  }
  clear_meshpoints (&bestMeshPoints);

  if (newLeftSaddle != NULL) {
    if (newLeftSaddle->s)
      free (newLeftSaddle->s);
    free (newLeftSaddle);
  }
  if (newRightSaddle != NULL) {
    if (newRightSaddle->s)
      free (newRightSaddle->s);
    free (newRightSaddle);
  }
  free (Saddle->s);
  free (Saddle);
  free_path (foldingPath);
  freeRNAWalkArrays ();
}

/**
 * Compute a distance based path between structure s1 and s2.
 * @param seq - the RNA sequence (ACGU)
 * @param s1 - first structure in dot-bracket format.
 * @param s2 - second structure in dot-bracket format.
 * @param iteration - the current iteration in a series of recursive calls. It should be initialized with 0.
 * @param maxIterations - maximal iterations of structurewalks.
 * @param maxKeep - how many structures are being kept (get_path in findpath.c)
 * @param maxStorage - maximal number of meshpoints (see insert_meshpoint in meshpoint.c)
 */
vrna_path_t *
levelSaddlePoint2(const char *seq, const char *s1, const char *s2/*, int *num_entry*/, int iteration, int maxIterations, int maxKeep, int maxStorage,
                  int maximum_distance1, int maximum_distance2)
{

  int i;
  vrna_path_t *foldingPath = get_path (seq, s1, s2, maxKeep/*, &steps, circ*/);
  vrna_path_t *Saddle = getSaddlePoint (foldingPath/*, steps*/);

  float oldSaddleEn = Saddle->en;
  vrna_path_t *newLeftSaddle, *newRightSaddle;
  vrna_path_t *path_left = NULL;
  vrna_path_t *path_right = NULL;
  vrna_path_t *newPath = NULL;
  unsigned int steps1, steps2;
  float newSaddleEn = oldSaddleEn;
  meshpoint_list bestMeshPoints;
  init_meshpoint_list (&bestMeshPoints);

  /* begin the nice pathfinding routine */
  short *pt1, *pt2;
  pt1 = vrna_ptable (s1);
  pt2 = vrna_ptable (s2);
  int n = pt1[0];

  short *intersect = (short *) vrna_alloc (sizeof(short) * (n + 2));
  intersect[0] = n;
  /* compute symmetrical difference between both structures */
  int a = 0;
  int b = 0;
  for (i = 1; i <= n; i++) {
    if (pt1[i] == pt2[i]) {
      intersect[i] = pt1[i];
      pt1[i] = pt2[i] = 0;
    }
    else {
      if (i < pt1[i])
        a++;
      if (i < pt2[i])
        b++;
      intersect[i] = 0;
    }
  }

  /* first collect all meshpoints where we get an initially better path */
  if (a + b > 1) {
    vrna_md_t md;
    vrna_md_set_default (&md);
    md.circ = circ;
    md.uniq_ML = 1;

    vrna_fold_compound_t *vc = vrna_fold_compound_TwoD (seq, s1, s2, &md, VRNA_OPTION_MFE);
    vrna_sol_TwoD_t *mfe_s = vrna_mfe_TwoD (vc, a + maximum_distance1, b + maximum_distance2);

    vrna_fold_compound_free (vc);

    for (i = 0; mfe_s[i].k != INF; i++) {
      if (mfe_s[i].k == -1) {
        free (mfe_s[i].s);
        continue;
      }
      path_left = get_path (seq, s1, mfe_s[i].s, maxKeep/*, &steps1, circ*/);
      path_right = get_path (seq, mfe_s[i].s, s2, maxKeep/*, &steps2, circ*/);

      newLeftSaddle = getSaddlePoint (path_left/*, steps1*/);
      newRightSaddle = getSaddlePoint (path_right/*, steps2*/);
      newSaddleEn = MAX2(newLeftSaddle->en, newRightSaddle->en);
      if (newSaddleEn <= oldSaddleEn) {
        insert_meshpoint_with_struct_energy (mfe_s[i].s, newSaddleEn, mfe_s[i].en, &bestMeshPoints, maxStorage);
      }
      free (mfe_s[i].s);
      free_path (path_left);
      path_left = NULL;
      free_path (path_right);
      path_right = NULL;
    }
    free (mfe_s);

#if 0
    for(i = 0; i<= a + maximum_distance1; i++) {
      for(j = 0; j<= b + maximum_distance2; j++) {
        if(dfold_structs[i][j].en != (float)INF/100.) {
          path_left = get_path(seq, s1, dfold_structs[i][j].s, maxKeep/*, &steps1, circ*/);
          path_right = get_path(seq, dfold_structs[i][j].s, s2, maxKeep/*, &steps2, circ*/);

          newLeftSaddle = getSaddlePoint(path_left/*, steps1*/);
          newRightSaddle = getSaddlePoint(path_right/*, steps2*/);
          newSaddleEn = MAX2(newLeftSaddle->en, newRightSaddle->en);
          if(newSaddleEn <= oldSaddleEn) {
            insert_meshpoint_with_struct_energy(dfold_structs[i][j].s, newSaddleEn, dfold_structs[i][j].en, &bestMeshPoints, maxStorage);
          }
          free_path(path_left); path_left=NULL;
          free_path(path_right); path_right=NULL;
        }
      }
      free(dfold_structs[i]);
    }
#endif
  }
  if (bestMeshPoints.count > 0) {
    /*
     now as we know n better SaddlePoints, we can iterate deeper to maybe obtain an even better saddle
     */

    if (iteration < maxIterations) {
      meshpoint *cur;
      /* int t_steps1, t_steps2; */
      vrna_path_t *t_path_left = NULL, *t_path_right = NULL;
      path_left = NULL;
      path_right = NULL;
      for (cur = bestMeshPoints.first; cur != NULL; cur = cur->next) {
        t_path_left = levelSaddlePoint2 (seq, s1, cur->s, /*&t_steps1,*/iteration + 1, maxIterations, maxKeep, maxStorage, maximum_distance1,
                                         maximum_distance2);
        newLeftSaddle = getSaddlePoint (t_path_left/*, t_steps1*/);

        t_path_right = levelSaddlePoint2 (seq, cur->s, s2, /*&t_steps2,*/iteration + 1, maxIterations, maxKeep, maxStorage, maximum_distance1,
                                          maximum_distance2);
        newRightSaddle = getSaddlePoint (t_path_right/*, t_steps2*/);

        newSaddleEn = MAX2(newLeftSaddle->en, newRightSaddle->en);

        if (newSaddleEn < oldSaddleEn) {
          if (path_left)
            free_path (path_left);
          path_left = t_path_left;
          if (path_right)
            free_path (path_right);
          path_right = t_path_right;
          /* steps1 = t_steps1; */
          /* steps2 = t_steps2; */
          oldSaddleEn = newSaddleEn;
        }
        else {
          free_path (t_path_left);
          t_path_left = NULL;
          free_path (t_path_right);
          t_path_right = NULL;
        }
      }
      if (path_left == NULL || path_right == NULL) {
        path_left = get_path (seq, s1, (bestMeshPoints.first)->s, maxKeep/*, &steps1, circ*/);
        path_right = get_path (seq, (bestMeshPoints.first)->s, s2, maxKeep/*, &steps2, circ*/);

      }
    }
    /* if we are in the last iteration step, we just take the best found in this round... */
    else {
      path_left = get_path (seq, s1, (bestMeshPoints.first)->s, maxKeep/*, &steps1, circ*/);
      path_right = get_path (seq, (bestMeshPoints.first)->s, s2, maxKeep/*, &steps2, circ*/);

      newLeftSaddle = getSaddlePoint (path_left/*, steps1*/);
      newRightSaddle = getSaddlePoint (path_right/*, steps2*/);
      newSaddleEn = MAX2(newLeftSaddle->en, newRightSaddle->en);

    }
    clear_meshpoints (&bestMeshPoints);
    for (steps1 = 0; path_left[steps1].s; steps1++)
      ;
    for (steps2 = 0; path_right[steps2].s; steps2++)
      ;
    newPath = (vrna_path_t *) vrna_alloc ((steps1 + steps2) * sizeof(vrna_path_t));
    memcpy ((vrna_path_t *) newPath, (vrna_path_t *) path_left, steps1 * sizeof(vrna_path_t));
    if (steps2 > 0) {
      memcpy (((vrna_path_t *) newPath) + (steps1), ((vrna_path_t *) path_right) + 1, (steps2 - 1) * sizeof(vrna_path_t));
      free (path_right[0].s); /* since we skipped this entry and it never would be free'd */
    }
  }
  else {
    free (pt1);
    free (pt2);
    free (intersect);
    return foldingPath;
  }

  free (pt1);
  free (pt2);
  free (intersect);
  free_path (foldingPath);
  if (path_left)
    free (path_left); /* do not free the structures, since they are further uses in newPath */
  if (path_right)
    free (path_right); /* do not free the structures, since they are further uses in newPath */

  //fprintf(stdout, "\rsearching for alternative paths...(%2.1f %% done %d %d)", (float)(curr_iteration)/(float)pow(2*maxStorage, maxIterations) * 100.,curr_iteration, (int)pow(2*maxStorage, maxIterations)  );
  fprintf (stdout, ".");
  fflush (stdout);

  return newPath;
}


static unsigned int
find_saddle_point(vrna_path_t *folding_path)
{
  unsigned int  i, position = 0;
  double        max_e, curr_e;

  max_e = folding_path[0].en;

  for (i = 1; folding_path[i].s; i++) {
    curr_e = folding_path[i].en;
    if (curr_e > max_e) {
      max_e     = curr_e;
      position  = i;
    }
  }

  return position;
}


static vrna_path_t *
clone_path_element(vrna_path_t *source)
{
  vrna_path_t *copy = vrna_alloc(2 * sizeof(vrna_path_t));

  copy[0].s   = strdup(source->s);
  copy[0].en  = source->en;
  copy[1].s   = NULL;
  copy[1].en  = 0;

  return copy;
}


static vrna_path_t *
concat_path(vrna_path_t *target,
            vrna_path_t *extension)
{
  /* find out size of both paths */
  unsigned int l1, l2;
  vrna_path_t  *ptr;

  for (l1 = 0, ptr = target; ptr->s; ptr++, l1++);
  for (l2 = 0, ptr = extension; ptr->s; ptr++, l2++);

  /* resize target path */
  target = (vrna_path_t *)vrna_realloc(target, sizeof(vrna_path_t) * (l1 + l2 + 1));

  memcpy(target + l1, extension, l2 * sizeof(vrna_path_t));

  target[l1 + l2].s   = NULL;
  target[l1 + l2].en  = 0;

  return target;
}


vrna_path_t *
getSaddlePoint(vrna_path_t *foldingPath)
{
  if (foldingPath)
    return clone_path_element(foldingPath + find_saddle_point(foldingPath));

  return NULL;
}


typedef struct {
  double  saddle_en;
  char    *saddle_structure;
} heap_elem_t;

typedef struct {
  unsigned int  num_elem;
  unsigned int  storage_size;
  heap_elem_t   **data;
} min_heap_t;


static min_heap_t *
init_heap(unsigned int max_nodes)
{
  min_heap_t  *heap = (min_heap_t *)vrna_alloc(sizeof(min_heap_t));

  heap->num_elem      = 0;
  heap->storage_size  = max_nodes;
  heap->data          = (heap_elem_t **)vrna_alloc(sizeof(heap_elem_t *) * (heap->storage_size + 1));

  return heap;
}

static void
destroy_heap(min_heap_t *heap)
{
  free(heap->data);
  free(heap);
}


void heap_insert(min_heap_t  *heap,
            heap_elem_t *d)
{
  if (heap->num_elem < heap->storage_size) {
    heap->data[heap->num_elem] = d;

    /* restore min-heap property */
    if (heap->num_elem > 0) {
      unsigned int current = heap->num_elem;
      unsigned int parent = floor((current - 1) / 2.);
      while (parent >= 0) {
        if (heap->data[parent]->saddle_en <= heap->data[current]->saddle_en)
          break;

        /* swap data with parent */
        heap_elem_t *tmp = heap->data[parent];
        heap->data[parent] = heap->data[current];
        heap->data[heap->num_elem] = tmp;

        current = parent;
        parent  = floor((current - 1) / 2.);
      }
    }
  }

}
