/* gcc -fopenmp -g3 -DTEST_FINDPATH findpath.c -o FINDpath -lRNA -lm -I ../H/ -L ./ */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/findpath.h"
#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/cofold.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/structure_utils.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define LOOP_EN

/**
 *  @brief
 */
typedef struct move {
  int i;  /* i,j>0 insert; i,j<0 delete */
  int j;
  int when;  /* 0 if still available, else resulting distance from start */
  int E;
} move_t;

/**
 *  @brief
 */
typedef struct intermediate {
  short *pt;      /**<  @brief  pair table */
  int Sen;        /**<  @brief  saddle energy so far */
  int curr_en;    /**<  @brief  current energy */
  move_t *moves;  /**<  @brief  remaining moves to target */
} intermediate_t;

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/

/*
#################################
# PRIVATE VARIABLES             #
#################################
*/
PRIVATE int         BP_dist;
PRIVATE move_t      *path=NULL;
PRIVATE int         path_fwd; /* 1: struc1->struc2, else struc2 -> struc1 */

PRIVATE vrna_fold_compound_t  *backward_compat_compound = NULL;

#ifdef _OPENMP

/* NOTE: all variables are assumed to be uninitialized if they are declared as threadprivate
*/
#pragma omp threadprivate(BP_dist, path, path_fwd, backward_compat_compound)

#endif

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE move_t  *copy_moves(move_t *mvs);
PRIVATE int     compare_ptable(const void *A, const void *B);
PRIVATE int     compare_energy(const void *A, const void *B);
PRIVATE int     compare_moves_when(const void *A, const void *B);
PRIVATE void    free_intermediate(intermediate_t *i);

#ifdef TEST_FINDPATH

/* TEST_FINDPATH, COFOLD */
PRIVATE void  usage(void);

#endif

PRIVATE int     find_path_once(vrna_fold_compound_t *vc, const char *struc1, const char *struc2, int maxE, int maxl);
PRIVATE int     try_moves(vrna_fold_compound_t *vc, intermediate_t c, int maxE, intermediate_t *next, int dist);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/
PUBLIC void free_path(vrna_path_t *path){
  vrna_path_t *tmp = path;
  if(tmp){
    while(tmp->s){ free(tmp->s); tmp++;}
    free(path);
  }
}

PUBLIC int
find_saddle(const char *seq,
            const char *struc1,
            const char *struc2,
            int max){

  int                 maxE;
  char                *sequence;
  vrna_fold_compound_t  *vc;
  vrna_md_t           md;

  vc = NULL;
  set_model_details(&md);

  if(backward_compat_compound){
    if(!strcmp(seq, backward_compat_compound->sequence)){ /* check if sequence is the same as before */
      md.window_size = backward_compat_compound->length;
      if(!memcmp(&md, &(backward_compat_compound->params->model_details), sizeof(vrna_md_t))){ /* check if model_details are the same as before */
        vc = backward_compat_compound; /* re-use previous vrna_fold_compound_t */
      }
    }
  }

  if(!vc){
    vrna_fold_compound_free(backward_compat_compound);

    sequence = vrna_cut_point_insert(seq, cut_point);

    backward_compat_compound = vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);

    free(sequence);
  }

  maxE = vrna_path_findpath_saddle(vc, struc1, struc2, max);

  return maxE;
}

PUBLIC int
vrna_path_findpath_saddle( vrna_fold_compound_t *vc,
                  const char *struc1,
                  const char *struc2,
                  int max) {

  int         maxl, maxE;
  const char  *tmp;
  move_t      *bestpath=NULL;
  int         dir;

  path_fwd  = 0;
  maxE      = INT_MAX - 1;

  maxl = 1;
  do {
    int saddleE;
    path_fwd = !path_fwd;
    if (maxl>max) maxl=max;
    if(path) free(path);
    saddleE  = find_path_once(vc, struc1, struc2, maxE, maxl);
    if (saddleE<maxE) {
      maxE = saddleE;
      if (bestpath) free(bestpath);
      bestpath = path;
      path = NULL;
      dir = path_fwd;
    } else{
      free(path);path=NULL;
    }
    tmp=struc1;
    struc1=struc2;
    struc2=tmp;
    maxl *=2;
  } while (maxl<2*max);

  /* (re)set some globals */
  path=bestpath;
  path_fwd = dir;

  return maxE;
}

PUBLIC vrna_path_t *
get_path( const char *seq,
          const char *s1,
          const char* s2,
          int maxkeep){

  vrna_path_t              *route    = NULL;
  char                *sequence = NULL;
  vrna_fold_compound_t  *vc       = NULL;
  vrna_md_t           md;

  set_model_details(&md);

  if(backward_compat_compound){
    if(!strcmp(seq, backward_compat_compound->sequence)){ /* check if sequence is the same as before */
      if(!memcmp(&md, &(backward_compat_compound->params->model_details), sizeof(vrna_md_t))){ /* check if model_details are the same as before */
        vc = backward_compat_compound; /* re-use previous vrna_fold_compound_t */
      }
    }
  }

  if(!vc){
    vrna_fold_compound_free(backward_compat_compound);

    sequence = vrna_cut_point_insert(seq, cut_point);

    backward_compat_compound = vc = vrna_fold_compound(sequence, &md, VRNA_OPTION_EVAL_ONLY);

    free(sequence);
  }

  route = vrna_path_findpath(vc, s1, s2, maxkeep);

  return (route);
}

PUBLIC vrna_path_t *
vrna_path_findpath(vrna_fold_compound_t *vc,
              const char *s1,
              const char* s2,
              int maxkeep){

  int E, d;
  vrna_path_t *route=NULL;

  E = vrna_path_findpath_saddle(vc, s1, s2, maxkeep);

  route = (vrna_path_t *)vrna_alloc((BP_dist+2)*sizeof(vrna_path_t));

  qsort(path, BP_dist, sizeof(move_t), compare_moves_when);

  if (path_fwd) {
    /* memorize start of path */
    route[0].s  = strdup(s1);
    route[0].en = vrna_eval_structure(vc, s1);

    for (d=0; d<BP_dist; d++) {
      int i,j;
      route[d+1].s = strdup(route[d].s);
      i = path[d].i; j=path[d].j;
      if (i<0) { /* delete */
        route[d+1].s[(-i)-1] = route[d+1].s[(-j)-1] = '.';
      } else {
        route[d+1].s[i-1] = '('; route[d+1].s[j-1] = ')';
      }
      route[d+1].en = path[d].E/100.0;
    }
  }
  else {
    /* memorize start of path */

    route[BP_dist].s  = strdup(s2);
    route[BP_dist].en = vrna_eval_structure(vc, s2);

    for (d=0; d<BP_dist; d++) {
      int i,j;
      route[BP_dist-d-1].s = strdup(route[BP_dist-d].s);
      i = path[d].i;
      j = path[d].j;
      if (i<0) { /* delete */
        route[BP_dist-d-1].s[(-i)-1] = route[BP_dist-d-1].s[(-j)-1] = '.';
      } else {
        route[BP_dist-d-1].s[i-1] = '('; route[BP_dist-d-1].s[j-1] = ')';
      }
      route[BP_dist-d-1].en = path[d].E/100.0;
    }
  }

#if _DEBUG_FINDPATH_
  fprintf(stderr, "\n%s\n%s\n%s\n\n", seq, s1, s2);
  for (d=0; d<=BP_dist; d++)
    fprintf(stderr, "%s %6.2f\n", route[d].s, route[d].en);
  fprintf(stderr, "%d\n", *num_entry);
#endif

  free(path);path=NULL;
  return (route);
}

PRIVATE int
try_moves(vrna_fold_compound_t *vc,
          intermediate_t c,
          int maxE,
          intermediate_t *next,
          int dist){

  int *loopidx, len, num_next=0, en, oldE;
  move_t *mv;
  short *pt;

  len = c.pt[0];
  loopidx = vrna_loopidx_from_ptable(c.pt);
  oldE = c.Sen;
  for (mv=c.moves; mv->i!=0; mv++) {
    int i,j;
    if (mv->when>0) continue;
    i = mv->i; j = mv->j;
    pt = (short *) vrna_alloc(sizeof(short)*(len+1));
    memcpy(pt, c.pt,(len+1)*sizeof(short));
    if (j<0) { /*it's a delete move */
      pt[-i]=0;
      pt[-j]=0;
    } else { /* insert move */
      if ((loopidx[i] == loopidx[j]) && /* i and j belong to same loop */
          (pt[i] == 0) && (pt[j]==0)     /* ... and are unpaired */
          ) {
        pt[i]=j;
        pt[j]=i;
      } else {
        free(pt);
        continue; /* llegal move, try next; */
      }
    }
#ifdef LOOP_EN
    en = c.curr_en + vrna_eval_move_pt(vc, c.pt, i, j);
#else
    en = vrna_eval_structure_pt(vc, pt);
#endif
    if (en<maxE) {
      next[num_next].Sen = (en>oldE)?en:oldE;
      next[num_next].curr_en = en;
      next[num_next].pt = pt;
      mv->when=dist;
      mv->E = en;
      next[num_next++].moves = copy_moves(c.moves);
      mv->when=0;
    }
    else free(pt);
  }
  free(loopidx);
  return num_next;
}

PRIVATE int find_path_once(vrna_fold_compound_t *vc, const char *struc1, const char *struc2, int maxE, int maxl) {
  short *pt1, *pt2;
  move_t *mlist;
  int i, len, d, dist=0, result;
  intermediate_t *current, *next;

  pt1 = vrna_ptable(struc1);
  pt2 = vrna_ptable(struc2);
  len = (int) strlen(struc1);

  mlist = (move_t *) vrna_alloc(sizeof(move_t)*len); /* bp_dist < n */

  for (i=1; i<=len; i++) {
    if (pt1[i] != pt2[i]) {
      if (i<pt1[i]) { /* need to delete this pair */
        mlist[dist].i = -i;
        mlist[dist].j = -pt1[i];
        mlist[dist++].when = 0;
      }
      if (i<pt2[i]) { /* need to insert this pair */
        mlist[dist].i = i;
        mlist[dist].j = pt2[i];
        mlist[dist++].when = 0;
      }
    }
  }
  free(pt2);
  BP_dist = dist;
  current = (intermediate_t *) vrna_alloc(sizeof(intermediate_t)*(maxl+1));
  current[0].pt = pt1;
  current[0].Sen = current[0].curr_en = vrna_eval_structure_pt(vc, pt1);
  current[0].moves = mlist;
  next = (intermediate_t *) vrna_alloc(sizeof(intermediate_t)*(dist*maxl+1));

  for (d=1; d<=dist; d++) { /* go through the distance classes */
    int c, u, num_next=0;
    intermediate_t *cc;

    for (c=0; current[c].pt != NULL; c++) {
      num_next += try_moves(vc, current[c], maxE, next+num_next, d);
    }
    if (num_next==0) {
      for (cc=current; cc->pt != NULL; cc++) free_intermediate(cc);
      current[0].Sen=INT_MAX;
      break;
    }
    /* remove duplicates via sort|uniq
       if this becomes a bottleneck we can use a hash instead */
    qsort(next, num_next, sizeof(intermediate_t),compare_ptable);
    for (u=0,c=1; c<num_next; c++) {
      if (memcmp(next[u].pt,next[c].pt,sizeof(short)*len)!=0) {
        next[++u] = next[c];
      } else {
        free_intermediate(next+c);
      }
    }
    num_next = u+1;
    qsort(next, num_next, sizeof(intermediate_t),compare_energy);
    /* free the old stuff */
    for (cc=current; cc->pt != NULL; cc++) free_intermediate(cc);
    for (u=0; u<maxl && u<num_next; u++) {
      current[u] = next[u];
    }
    for (; u<num_next; u++)
      free_intermediate(next+u);
    num_next=0;
  }
  free(next);
  path = current[0].moves;
  result = current[0].Sen;
  free(current[0].pt); free(current);
  return(result);
}

PRIVATE void free_intermediate(intermediate_t *i) {
   free(i->pt);
   free(i->moves);
   i->pt = NULL;
   i->moves = NULL;
   i->Sen = INT_MAX;
 }

PRIVATE int compare_ptable(const void *A, const void *B) {
  intermediate_t *a, *b;
  int c;
  a = (intermediate_t *) A;
  b = (intermediate_t *) B;

  c = memcmp(a->pt, b->pt, a->pt[0]*sizeof(short));
  if (c!=0) return c;
  if ((a->Sen - b->Sen) != 0) return (a->Sen - b->Sen);
  return (a->curr_en - b->curr_en);
}

PRIVATE int compare_energy(const void *A, const void *B) {
  intermediate_t *a, *b;
  a = (intermediate_t *) A;
  b = (intermediate_t *) B;

  if ((a->Sen - b->Sen) != 0) return (a->Sen - b->Sen);
  return (a->curr_en - b->curr_en);
}

PRIVATE int compare_moves_when(const void *A, const void *B) {
  move_t *a, *b;
  a = (move_t *) A;
  b = (move_t *) B;

  return(a->when - b->when);
}

PRIVATE move_t* copy_moves(move_t *mvs) {
  move_t *new;
  new = (move_t *) vrna_alloc(sizeof(move_t)*(BP_dist+1));
  memcpy(new,mvs,sizeof(move_t)*(BP_dist+1));
  return new;
}

#ifdef TEST_FINDPATH

PUBLIC void print_path(const char *seq, const char *struc) {
  int d;
  char *s;
  s = strdup(struc);
  if (cut_point == -1)
    printf("%s\n%s\n", seq, s);
    /* printf("%s\n%s %6.2f\n", seq, s, vrna_eval_structure_simple(seq,s)); */
  else {
    char *pstruct, *pseq;
    pstruct = vrna_cut_point_insert(s, cut_point);
    pseq = vrna_cut_point_insert(seq, cut_point);
    printf("%s\n%s\n", pseq, pstruct);
    /* printf("%s\n%s %6.2f\n", pseq, pstruct, vrna_eval_structure_simple(seq,s)); */
    free(pstruct);
    free(pseq);
  }
  qsort(path, BP_dist, sizeof(move_t), compare_moves_when);
  for (d=0; d<BP_dist; d++) {
    int i,j;
    i = path[d].i; j=path[d].j;
    if (i<0) { /* delete */
      s[(-i)-1] = s[(-j)-1] = '.';
    } else {
      s[i-1] = '('; s[j-1] = ')';
    }
    /* printf("%s %6.2f - %6.2f\n", s, vrna_eval_structure_simple(seq,s), path[d].E/100.0); */
  }
  free(s);
}

int main(int argc, char *argv[]) {
  char *line, *seq, *s1, *s2;
  int E, maxkeep=1000;
  int verbose=0, i;
  vrna_path_t *route, *r;

  for (i=1; i<argc; i++) {
    switch ( argv[i][1] ) {
      case 'm': if (strcmp(argv[i],"-m")==0)
                  sscanf(argv[++i], "%d", &maxkeep);
                break;
      case 'v': verbose = !strcmp(argv[i],"-v");
                break;
      case 'd': if (strcmp(argv[i],"-d")==0){
                  sscanf(argv[++i], "%d", &dangles);
                  md.dangles = dangles;
                }
                break;
      default: usage();
    }
  }

  cut_point = -1;
  line = vrna_read_line(stdin);
  seq = vrna_cut_point_remove(line, &cut_point);
  free(line);   
  line = vrna_read_line(stdin);
  s1 = vrna_cut_point_remove(line, &cut_point);
  free(line);
  line = vrna_read_line(stdin);
  s2 = vrna_cut_point_remove(line, &cut_point);
  free(line);

  E = find_saddle(seq, s1, s2, maxkeep);
  printf("saddle_energy = %6.2f\n", E/100.);
  if (verbose) {
      if (path_fwd)
          print_path(seq,s1);
      else
          print_path(seq,s2);
      free(path);
      path = NULL;
      route = get_path(seq, s1, s2, maxkeep);
      for (r=route; r->s; r++) {
          if (cut_point == -1) {
              printf("%s %6.2f\n", r->s, r->en);
              /* printf("%s %6.2f - %6.2f\n", r->s, vrna_eval_structure_simple(seq,r->s), r->en); */
          } else {
              char *pstruct;
              pstruct = vrna_cut_point_insert(r->s, cut_point);
              printf("%s %6.2f\n", pstruct, r->en);
              /* printf("%s %6.2f - %6.2f\n", pstruct, vrna_eval_structure_simple(seq,r->s), r->en); */
              free(pstruct);
          }
          free(r->s);
      }
      free(route);
  }
  free(seq); free(s1); free(s2);
  return(EXIT_SUCCESS);
}

static void usage(void){
  vrna_message_error("usage: findpath.c  [-m depth] [-d[0|1|2]] [-v]");
}

#endif
