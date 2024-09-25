/*
        PostScript and GML output for RNA secondary structures
                    and pair probability matrices

                 c  Ivo Hofacker and Peter F Stadler
                          Vienna RNA package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "ViennaRNA/model.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/sequences/alignments.h"
#include "ViennaRNA/utils/log.h"
#include "ViennaRNA/plotting/probabilities.h"

#include "ViennaRNA/static/templates_postscript.h"

#include "ViennaRNA/plotting/ps_helpers.inc"


#ifndef INLINE
#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif
#endif

/*
#################################
# PRIVATE MACROS                #
#################################
*/

#define SIZE 452.
#define PMIN 0.00001

typedef struct {
  vrna_data_lin_t **data;
  char            **titles;
  size_t          size;
  size_t          mem;
} lin_data_container;

typedef struct {
  lin_data_container north;
  lin_data_container east;
  lin_data_container south;
  lin_data_container west;
} lin_dat;

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
#define dp_add_lindata(data_var, data, title_var, title, size, avail)  \
  (title_var)[size] = title; \
  (data_var)[size]  = data; \
  if((++size) == (avail)){ \
    avail *= 1.2; \
    data_var  = (vrna_data_lin_t **)vrna_realloc(data_var, sizeof(vrna_data_lin_t *) * avail); \
    title_var = (char **)vrna_realloc(title_var, sizeof(char *) * avail); \
  }


#define dp_finalize_lindata(data_var, title_var, size) \
  (data_var)[size]  = NULL; \
  (title_var)[size] = NULL; \
  (data_var)        = (vrna_data_lin_t **)vrna_realloc(data_var, sizeof(vrna_data_lin_t *) * (size + 1)); \
  (title_var)       = (char **)vrna_realloc(title_var, sizeof(char *) * (size + 1));

PRIVATE char *comment_dotplot = "This file contains the square roots of probabilities in the form\n"
                                "i  j  sqrt(p(i,j)) ubox";


/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE FILE  *PS_dot_common(const char *seq, unsigned int *nicks, const char *wastlfile, char *comment, int winsize, unsigned int options);
PRIVATE int   sort_plist_by_type_desc(const void *p1, const void *p2);
PRIVATE int   sort_plist_by_prob_asc(const void *p1, const void *p2);
PRIVATE int   sort_cpair_by_type_desc(const void *p1, const void *p2);
PRIVATE int   sort_cpair_by_prob_asc(const void *p1, const void *p2);
PRIVATE void  EPS_print_title(FILE *eps, const char *title);
PRIVATE void  EPS_print_header(FILE *eps, int bbox[4], const char *comment, unsigned int options);
PRIVATE void  EPS_print_ud_data(FILE *eps, plist *pl, plist *mf);
PRIVATE void  EPS_print_sc_motif_data(FILE *eps, plist *pl, plist *mf);
PRIVATE void  EPS_print_bpp_data(FILE *eps, plist *pl, plist *mf);
PRIVATE void  EPS_print_linear_data(FILE *eps, const char *varname, lin_data_container *data);
PRIVATE vrna_data_lin_t *plist_to_accessibility(plist *pl, unsigned int length);
PRIVATE vrna_data_lin_t *plist_to_ud_motif_prob(plist *pl, unsigned int length);
PRIVATE void  EPS_print_sd_data(FILE *eps, vrna_ep_t *pl, vrna_ep_t *mf);

INLINE PRIVATE int
init_lin_data(lin_dat *dat);

INLINE PRIVATE void
free_lin_data(lin_dat *dat);

INLINE PRIVATE int
push_lin_data(lin_data_container  *c,
              vrna_data_lin_t     *data,
              char                *title);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC int
PS_color_dot_plot(char *seq,
                  cpair *pi,
                  char *wastlfile){

  /* produce color PostScript dot plot from cpair */

  FILE          *wastl;
  unsigned int  *nicks;
  int           i, gq_num, pi_size;
  cpair         *ptr;

  nicks = NULL;

  if (cut_point > 0) {
    nicks = (unsigned int *)vrna_alloc(sizeof(unsigned int) * 2);
    nicks[0] = cut_point;
    nicks[1] = 0;
  }

  wastl = PS_dot_common(seq, nicks, wastlfile, NULL, 0, PS_MACRO_DOTPLOT_SD);

  free(nicks);

  if (wastl==NULL)
    return 0; /* return 0 for failure */

  fprintf(wastl, "/hsb {\n"
          "dup 0.3 mul 1 exch sub sethsbcolor\n"
          "} bind def\n\n");

  fprintf(wastl,  "\n%%draw the grid\ndrawgrid\n\n");
  fprintf(wastl,"%%start of base pair probability data\n");

  for(gq_num = pi_size = 0, ptr = pi; ptr->i > 0; ptr++, pi_size++)
    if (ptr->type == VRNA_PLIST_TYPE_GQUAD) gq_num++;
  qsort(pi, pi_size, sizeof(cpair), sort_cpair_by_type_desc);
  /* sort all gquad triangles by probability to bring lower probs to the front */
  qsort(pi, gq_num, sizeof(cpair), sort_cpair_by_prob_asc);

  /* print boxes */
  i=0;
  while (pi[i].j>0) {
    if(pi[i].type == VRNA_PLIST_TYPE_GQUAD){
      fprintf(wastl,
              "%d %d %1.6f utri\n",
              pi[i].i,
              pi[i].j,
              sqrt(pi[i].p));
    } else if ((pi[i].type == VRNA_PLIST_TYPE_BASEPAIR) ||
               (pi[i].type == VRNA_PLIST_TYPE_TRIPLE)) {
      fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.6f ubox\n",
             pi[i].hue,
             pi[i].sat,
             (pi[i].i < pi[i].j) ? pi[i].i : pi[i].j,
             (pi[i].i < pi[i].j) ? pi[i].j : pi[i].i,
             sqrt(pi[i].p));

      if (pi[i].mfe)
        fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.4f lbox\n",
               pi[i].hue,
               pi[i].sat,
               (pi[i].i < pi[i].j) ? pi[i].i : pi[i].j,
               (pi[i].i < pi[i].j) ? pi[i].j : pi[i].i,
               pi[i].p);
    }
    i++;
  }

  print_PS_footer(wastl);

  fclose(wastl);
  return 1; /* success */
}


PUBLIC int
PS_dot_plot_list( char *seq,
                  char *wastlfile,
                  plist *pl,
                  plist *mf,
                  char *comment){

  return vrna_plot_dp_PS_list(seq, cut_point, wastlfile, pl, mf, comment);
}


PUBLIC int
vrna_plot_dp_PS_list( char *seq,
                      int cp,
                      char *wastlfile,
                      plist *pl,
                      plist *mf,
                      char *comment){

  FILE          *wastl;
  size_t        cnt;
  char          *seq_plain, **seqs;
  unsigned int  *nicks, curr_length;
  int           pl_size, gq_num;
  plist         *pl1;

  nicks     = NULL;
  seq_plain = NULL;
  seqs      = vrna_strsplit(seq, "&");

  if (seqs) {
    /* count total number of strands */
    for (cnt = 0; seqs[cnt]; cnt++);

    seq_plain   = seqs[0];
    curr_length = strlen(seq_plain);

    if (seqs[1]) {
      nicks = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (cnt + 1));

      /* get first cut point and concatenate sequences */
      nicks[0] = curr_length + 1;
      vrna_strcat_printf(&seq_plain, "%s", seqs[1]);
      curr_length += strlen(seqs[1]);
      free(seqs[1]);

      /* add all remaining sequences (so far, we do not store the individual nicks for the remaining sequences */
      cnt = 2;
      while(seqs[cnt]) {
        nicks[cnt - 1] = curr_length + 1;
        vrna_strcat_printf(&seq_plain, "%s", seqs[cnt]);
        curr_length += strlen(seqs[cnt]);
        free(seqs[cnt++]);
      }
    }

    free(seqs);
  }

  wastl = PS_dot_common(seq_plain, nicks, wastlfile, comment, 0, PS_MACRO_DOTPLOT_ALL);

  free(seq_plain);
  free(nicks);

  if (wastl==NULL)
    return 0; /* return 0 for failure */

  fprintf(wastl,"%%data starts here\n");

  /* sort the plist to bring all gquad triangles to the front */
  if (pl) {
    for(gq_num = pl_size = 0, pl1 = pl; pl1->i > 0; pl1++, pl_size++)
      if (pl1->type == VRNA_PLIST_TYPE_GQUAD) gq_num++;
    qsort(pl, pl_size, sizeof(plist), sort_plist_by_type_desc);
    /* sort all gquad triangles by probability to bring lower probs to the front */
    qsort(pl, gq_num, sizeof(plist), sort_plist_by_prob_asc);
  }

  EPS_print_sd_data(wastl, pl, mf);
  EPS_print_sc_motif_data(wastl, pl, mf);

  fprintf(wastl, "\n%%draw the grid\ndrawgrid\n\n");
  fprintf(wastl,"%%start of base pair probability data\n");

  EPS_print_bpp_data(wastl, pl, mf);

  print_PS_footer(wastl);

  fclose(wastl);
  return 1; /* success */
}


PUBLIC int
vrna_plot_dp_EPS( const char              *filename,
                  const char              *sequence,
                  vrna_ep_t            *upper,
                  vrna_ep_t            *lower,
                  vrna_dotplot_auxdata_t  *auxdata,
                  unsigned int            options){

  char            **lintoptitle,**linbottomtitle,**linlefttitle,**linrighttitle,
                  *c, *comment, *title;
  int             i, lintop_num, lintop_size, linbottom_num,
                  linbottom_size, linleft_num, linleft_size, linright_num,
                  linright_size, bbox[4];
  FILE            *fh;
  vrna_data_lin_t **lintop, **linbottom, **linleft, **linright, *ud_lin, *access;
  lin_dat         linear_data;
  
  fh = fopen(filename, "w");
  if(!fh){
    vrna_log_warning("can't open %s for dot plot", filename);
    return 0; /* return 0 for failure */
  }

  comment = title = NULL;

  init_lin_data(&linear_data);

  lintop_num      = 0;
  lintop_size     = 5;
  linbottom_num   = 0;
  linbottom_size  = 5;
  linleft_num     = 0;
  linleft_size    = 5;
  linright_num    = 0;
  linright_size   = 5;
  bbox[0]         = 0;
  bbox[1]         = 0;
  bbox[2]         = 700;
  bbox[3]         = 720;
  access          = NULL;
  ud_lin          = NULL;
  lintop          = (vrna_data_lin_t **)vrna_alloc(sizeof(vrna_data_lin_t *) * lintop_size);
  lintoptitle     = (char **)vrna_alloc(sizeof(char *) * lintop_size);
  linbottom       = (vrna_data_lin_t **)vrna_alloc(sizeof(vrna_data_lin_t *) * linbottom_size);
  linbottomtitle  = (char **)vrna_alloc(sizeof(char *) * linbottom_size);
  linleft         = (vrna_data_lin_t **)vrna_alloc(sizeof(vrna_data_lin_t *) * linleft_size);
  linlefttitle    = (char **)vrna_alloc(sizeof(char *) * linleft_size);
  linright        = (vrna_data_lin_t **)vrna_alloc(sizeof(vrna_data_lin_t *) * linright_size);
  linrighttitle   = (char **)vrna_alloc(sizeof(char *) * linright_size);

  /* prepare linear data and retrieve number of additional linear data lines to correct bounding box */
  if(options & VRNA_PLOT_PROBABILITIES_UD_LIN){
    ud_lin = plist_to_ud_motif_prob(upper, strlen(sequence));
    if(ud_lin){
      push_lin_data(&(linear_data.north), ud_lin, "Protein binding");
      push_lin_data(&(linear_data.east), ud_lin, "Protein binding");
      push_lin_data(&(linear_data.south), ud_lin, "Protein binding");
      push_lin_data(&(linear_data.west), ud_lin, "Protein binding");
    }
  }

  if(options & VRNA_PLOT_PROBABILITIES_ACC){
    access = plist_to_accessibility(upper, strlen(sequence));
    push_lin_data(&(linear_data.north), access, "Accessibility");
  }

  if(auxdata){
    if(auxdata->top){
      for(i = 0; auxdata->top[i]; i++){
        push_lin_data(&(linear_data.north), auxdata->top[i], auxdata->top_title[i]);
      }
    }
    if(auxdata->bottom){
      for(i = 0; auxdata->bottom[i]; i++){
        push_lin_data(&(linear_data.south), auxdata->bottom[i], auxdata->bottom_title[i]);
      }
    }
    if(auxdata->left){
      for(i = 0; auxdata->left[i]; i++){
        push_lin_data(&(linear_data.west), auxdata->left[i], auxdata->left_title[i]);
      }
    }
    if(auxdata->right){
      for(i = 0; auxdata->right[i]; i++){
        push_lin_data(&(linear_data.east), auxdata->right[i], auxdata->right_title[i]);
      }
    }
  }

  if(auxdata){
    comment = auxdata->comment;
    title   = (auxdata->title) ? strdup(auxdata->title) : NULL;
  }

  if(!title){
    title = strdup(filename);
    if((c=strrchr(title, '_')))
      *c='\0';
  }

  /* start printing postscript file */
  EPS_print_header(fh, bbox, comment, PS_MACRO_DOTPLOT_ALL);

  EPS_print_title(fh, title);
  print_PS_sequence(fh, sequence);

  fprintf(fh,"%% BEGIN linear data array\n\n");
  EPS_print_linear_data(fh, "topData", &(linear_data.north));
  EPS_print_linear_data(fh, "leftData", &(linear_data.west));
  EPS_print_linear_data(fh, "bottomData", &(linear_data.south));
  EPS_print_linear_data(fh, "rightData", &(linear_data.east));

  fprintf(fh,"%% END linear data arrays\n");

  fprintf(fh, "\n%%Finally, prepare canvas\n\n"
              "%%draw title\ndrawTitle\n\n"
              "%%prepare coordinate system, draw grid and sequence\n"
              "/Helvetica findfont 0.95 scalefont setfont\n\n"
              "%%prepare coordinate system\nprepareCoords\n\n"
              "%%draw sequence arround grid\ndrawseq\n\n"
              "%%draw grid\ndrawgrid\n\n"
              "%%draw auxiliary linear data (if available)\ndrawData\n\n");

  fprintf(fh,"%%data (commands) starts here\n");

  if(options & VRNA_PLOT_PROBABILITIES_SD){
    EPS_print_sd_data(fh, upper, lower);
  }

  if(options & VRNA_PLOT_PROBABILITIES_UD){
    EPS_print_ud_data(fh, upper, lower);
  }


  EPS_print_sc_motif_data(fh, upper, lower);
  EPS_print_bpp_data(fh, upper, lower);

  print_PS_footer(fh);

  free_lin_data(&linear_data);

  fclose(fh);
  free(lintoptitle);
  free(lintop);
  free(linbottomtitle);
  free(linbottom);
  free(linlefttitle);
  free(linleft);
  free(linrighttitle);
  free(linright);
  free(access);
  free(ud_lin);
  free(title);

  return 1; /* success */
}


PUBLIC int
PS_color_dot_plot_turn( char *seq,
                        cpair *pi,
                        char *wastlfile,
                        int winSize) {

  /* produce color PostScript dot plot from cpair */

  FILE          *wastl;
  unsigned int  *nicks;
  int           i;

  nicks = NULL;

  if (cut_point > 0) {
    nicks = (unsigned int *)vrna_alloc(sizeof(unsigned int) * 2);
    nicks[0] = cut_point;
    nicks[1] = 0;
  }

  wastl = PS_dot_common(seq, nicks, wastlfile, NULL, winSize, 0);

  free(nicks);

  if (wastl==NULL)
    return 0; /* return 0 for failure */

  fprintf(wastl, "/hsb {\n"
          "dup 0.3 mul 1 exch sub sethsbcolor\n"
          "} bind def\n\n"
          "%%BEGIN DATA\n");

  if(winSize > 0)
    fprintf(wastl, "\n%%draw the grid\ndrawgrid_turn\n\n");
  else
    fprintf(wastl,  "\n%%draw the grid\ndrawgrid\n\n");
  fprintf(wastl,"%%start of base pair probability data\n");

  /* print boxes */
   i=0;
   while (pi[i].j>0) {
     fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.6f ubox\n",
             pi[i].hue,
             pi[i].sat,
             (pi[i].i < pi[i].j) ? pi[i].i : pi[i].j,
             (pi[i].i < pi[i].j) ? pi[i].j : pi[i].i,
             sqrt(pi[i].p));

     if (pi[i].mfe)
       fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.4f lbox\n",
               pi[i].hue,
               pi[i].sat,
               (pi[i].i < pi[i].j) ? pi[i].i : pi[i].j,
               (pi[i].i < pi[i].j) ? pi[i].j : pi[i].i,
               pi[i].p);
     i++;
   }

   print_PS_footer(wastl);

   fclose(wastl);
   return 1; /* success */
}


PUBLIC int
PS_dot_plot_turn( char *seq,
                  plist *pl,
                  char *wastlfile,
                  int winSize) {

  /* produce color PostScript dot plot from cpair */

  FILE          *wastl;
  unsigned int  *nicks;
  int           i;

  nicks = NULL;

  if (cut_point > 0) {
    nicks = (unsigned int *)vrna_alloc(sizeof(unsigned int) * 2);
    nicks[0] = cut_point;
    nicks[1] = 0;
  }

  wastl = PS_dot_common(seq, nicks, wastlfile, NULL, winSize, 0);

  free(nicks);

  if (wastl==NULL)
    return 0; /* return 0 for failure */

  if(winSize > 0)
    fprintf(wastl, "\n%%draw the grid\ndrawgrid_turn\n\n");
  else
    fprintf(wastl,  "\n%%draw the grid\ndrawgrid\n\n");
  fprintf(wastl,"%%start of base pair probability data\n");
  /* print boxes */
  i=0;
  if (pl) {
    while (pl[i].j>0) {
      fprintf(wastl,"%d %d %1.4f ubox\n",
             (pl[i].i < pl[i].j) ? pl[i].i : pl[i].j,
             (pl[i].i < pl[i].j) ? pl[i].j : pl[i].i,
              sqrt(pl[i].p));
      i++;
    }
  }

  print_PS_footer(wastl);

  fclose(wastl);
  return 1; /* success */
}


/*
#####################################
# BEGIN OF STATIC HELPER FUNCTIONS  #
#####################################
*/

PRIVATE void
EPS_print_title(FILE *eps, const char *title){

  fprintf(eps,"/DPtitle {\n  (%s)\n} def\n\n", title);
}


PRIVATE void
EPS_print_header( FILE          *eps,
                  int           bbox[4],
                  const char    *comment,
                  unsigned int  options){

  char  *full_comment = NULL;

  vrna_md_t md;

  set_model_details(&md);

  if (comment)
    full_comment = vrna_strdup_printf("%s\n\n%s",
                                      comment,
                                      comment_dotplot);
  else
    full_comment = comment_dotplot;

  print_PS_header(eps,
                  "RNA Dot Plot",
                  bbox,
                  &md,
                  full_comment,
                  "DPdict",
                  PS_MACRO_DOTPLOT_BASE | options);

  if (comment)
    free(full_comment);
}


PRIVATE void
EPS_print_sd_data(FILE          *eps,
                  vrna_ep_t  *pl,
                  vrna_ep_t  *mf){

  int     pl_size, gq_num;
  double  tmp;
  plist   *pl1;

  /* sort the plist to bring all gquad triangles to the front */
  if (pl) {
    for(gq_num = pl_size = 0, pl1 = pl; pl1->i > 0; pl1++, pl_size++)
      if(pl1->type == VRNA_PLIST_TYPE_GQUAD) gq_num++;

    qsort(pl, pl_size, sizeof(plist), sort_plist_by_type_desc);

    /* sort all gquad triangles by probability to bring lower probs to the front */
    qsort(pl, gq_num, sizeof(plist), sort_plist_by_prob_asc);
  }

  /* print triangles for g-quadruplexes in upper half */
  fprintf(eps,"\n%%start of quadruplex data\n");

  if (pl)
    for (pl1=pl; pl1->i > 0; pl1++) {
      if(pl1->type == VRNA_PLIST_TYPE_GQUAD){
        tmp = sqrt(pl1->p);
        fprintf(eps, "%d %d %1.9f utri\n", pl1->i, pl1->j, tmp);
      }
    }
}


PRIVATE void
EPS_print_sc_motif_data(FILE          *eps,
                        vrna_ep_t  *pl,
                        vrna_ep_t  *mf){

  double  tmp;
  plist   *pl1;

  /* print triangles for hairpin loop motifs in upper half */
  fprintf(eps,"\n%%start of Hmotif data\n");
  if (pl)
    for (pl1=pl; pl1->i > 0; pl1++) {
      if(pl1->type == VRNA_PLIST_TYPE_H_MOTIF){
        tmp = sqrt(pl1->p);
        fprintf(eps, "%d %d %1.9f uHmotif\n", pl1->i, pl1->j, tmp);
      }
    }

  if (mf)
    for (pl1=mf; pl1->i > 0; pl1++) {
      if(pl1->type == VRNA_PLIST_TYPE_H_MOTIF){
        tmp = sqrt(pl1->p);
        fprintf(eps, "%d %d %1.9f lHmotif\n", pl1->i, pl1->j, tmp);
      }
    }

  /* print triangles for internal loop motifs in upper half */
  fprintf(eps,"\n%%start of Imotif data\n");
  int   a,b;
  float ppp;
  a = b = 0;

  if (pl)
    for (pl1=pl; pl1->i > 0; pl1++) {
      if(pl1->type == VRNA_PLIST_TYPE_I_MOTIF){
        if(a == 0){
          a = pl1->i;
          b = pl1->j;
          ppp = tmp = sqrt(pl1->p);
        } else {
          fprintf(eps, "%d %d %d %d %1.9f uImotif\n", a, b, pl1->i, pl1->j, ppp);
          a = b = 0;
        }
      }
    }

  if (mf)
    for (a = b= 0, pl1=mf; pl1->i > 0; pl1++) {
      if(pl1->type == VRNA_PLIST_TYPE_I_MOTIF){
        if(a == 0){
          a = pl1->i;
          b = pl1->j;
          ppp = sqrt(pl1->p);
        } else {
          fprintf(eps, "%d %d %d %d %1.9f lImotif\n", a, b, pl1->i, pl1->j, ppp);
          a = b = 0;
        }
      }
    }
}


PRIVATE void
EPS_print_bpp_data( FILE          *eps,
                    vrna_ep_t  *pl,
                    vrna_ep_t  *mf){

  double  tmp;
  plist   *pl1;

  fprintf(eps,"%%start of base pair probability data\n");

  /* print boxes in upper right half*/
  if (pl)
    for (pl1 = pl; pl1->i>0; pl1++) {
      tmp = sqrt(pl1->p);
      if ((pl1->type == VRNA_PLIST_TYPE_BASEPAIR) ||
          (pl1->type == VRNA_PLIST_TYPE_TRIPLE))
          fprintf(eps,
                  "%d %d %1.9f ubox\n",
                  (pl1->i < pl1->j) ? pl1->i : pl1->j,
                  (pl1->i < pl1->j) ? pl1->j : pl1->i,
                  tmp);
    }


  /* print boxes in lower left half (mfe) */
  if (mf)
    for (pl1=mf; pl1->i>0; pl1++) {
      tmp = sqrt(pl1->p);
      if ((pl1->type == VRNA_PLIST_TYPE_BASEPAIR) ||
          (pl1->type == VRNA_PLIST_TYPE_TRIPLE))
        fprintf(eps,
                "%d %d %1.7f lbox\n",
                (pl1->i < pl1->j) ? pl1->i : pl1->j,
                (pl1->i < pl1->j) ? pl1->j : pl1->i,
                tmp);
    }
}


PRIVATE void
EPS_print_ud_data(FILE          *eps,
                  vrna_ep_t  *pl,
                  vrna_ep_t  *mf){

  double  tmp;
  plist   *pl1;

  /* print triangles for unstructured domain motifs in upper half */
  fprintf(eps,"\n%%start of unstructured domain motif data\n");
  if (pl)
    for(pl1 = pl; pl1->i > 0; pl1++){
      if(pl1->type == VRNA_PLIST_TYPE_UD_MOTIF){
        tmp = sqrt(pl1->p);
        fprintf(eps, "%d %d %1.9f uUDmotif\n", pl1->i, pl1->j, tmp);
      }
    }

  if (mf)
    for(pl1 = mf; pl1->i > 0; pl1++){
      if(pl1->type == VRNA_PLIST_TYPE_UD_MOTIF){
        tmp = sqrt(pl1->p);
        fprintf(eps, "%d %d %1.9f lUDmotif\n", pl1->i, pl1->j, tmp);
      }
    }
}


PRIVATE void
EPS_print_linear_data(FILE                *eps,
                      const char          *varname,
                      lin_data_container  *c){

  size_t          i;
  vrna_data_lin_t *ptr;

  /* count number of data sets */

  fprintf(eps, "/%s [\n", varname);
  for(i = 0; i < c->size; i++){
    fprintf(eps, "[ (%s)\n", c->titles[i]);
    for(ptr = c->data[i]; ptr->position > 0; ptr++){
      if((ptr->color.hue + ptr->color.sat + ptr->color.bri) == 0.){
        fprintf(eps, "  [ %d %1.9f ]\n", ptr->position, ptr->value);
      } else {
        fprintf(eps, "  [ %d %1.9f %1.4f %1.4f %1.4f]\n", ptr->position, ptr->value, ptr->color.hue, ptr->color.sat, ptr->color.bri);
      }
    }
    fprintf(eps, "]\n");
  }
  fprintf(eps, "] def\n\n");
}


/*---------------------------------------------------------------------------*/


static int sort_plist_by_type_desc(const void *p1, const void *p2){
  if(((plist*)p1)->type > ((plist*)p2)->type) return -1;
  if(((plist*)p1)->type < ((plist*)p2)->type) return 1;

  /* same type?, order by position (ascending) */
  if (((plist*)p1)->i > ((plist*)p2)->i) return 1;
  if (((plist*)p1)->i < ((plist*)p2)->i) return -1;
  if (((plist*)p1)->j > ((plist*)p2)->j) return 1;
  if (((plist*)p1)->j < ((plist*)p2)->j) return -1;

  return 0;
}

static int sort_plist_by_prob_asc(const void *p1, const void *p2){
  if(((plist*)p1)->p > ((plist*)p2)->p) return 1;
  if(((plist*)p1)->p < ((plist*)p2)->p) return -1;

  /* same probability?, order by position (ascending) */
  if (((plist*)p1)->i > ((plist*)p2)->i) return 1;
  if (((plist*)p1)->i < ((plist*)p2)->i) return -1;
  if (((plist*)p1)->j > ((plist*)p2)->j) return 1;
  if (((plist*)p1)->j < ((plist*)p2)->j) return -1;
  return 0;
}


static int sort_cpair_by_type_desc(const void *p1, const void *p2){
  if(((cpair*)p1)->type > ((cpair*)p2)->type) return -1;
  if(((cpair*)p1)->type < ((cpair*)p2)->type) return 1;

  /* same type?, order by position (ascending) */
  if (((cpair*)p1)->i > ((cpair*)p2)->i) return 1;
  if (((cpair*)p1)->i < ((cpair*)p2)->i) return -1;
  if (((cpair*)p1)->j > ((cpair*)p2)->j) return 1;
  if (((cpair*)p1)->j < ((cpair*)p2)->j) return -1;

  return 0;
}

static int sort_cpair_by_prob_asc(const void *p1, const void *p2){
  if(((cpair*)p1)->p > ((cpair*)p2)->p) return 1;
  if(((cpair*)p1)->p < ((cpair*)p2)->p) return -1;

  /* same probability?, order by position (ascending) */
  if (((cpair*)p1)->i > ((cpair*)p2)->i) return 1;
  if (((cpair*)p1)->i < ((cpair*)p2)->i) return -1;
  if (((cpair*)p1)->j > ((cpair*)p2)->j) return 1;
  if (((cpair*)p1)->j < ((cpair*)p2)->j) return -1;
  return 0;
}


PRIVATE vrna_data_lin_t *
plist_to_accessibility(plist *pl, unsigned int length){

  int   i;
  plist *ptr;

  vrna_data_lin_t *data = NULL;

  data = (vrna_data_lin_t *)vrna_alloc(sizeof(vrna_data_lin_t) * (length + 1));

  for(ptr = pl; ptr->i > 0; ptr++){
    if(ptr->type == VRNA_PLIST_TYPE_BASEPAIR){
      data[ptr->i - 1].value += ptr->p;
      data[ptr->j - 1].value += ptr->p;
    }
  }

  /* go through the entire list and square-root probabilities again */
  for(i = 0; i < length; i++){
    data[i].position = i + 1;
    data[i].value    = sqrt((double)(1. - data[i].value));
  }

  data[length].position = 0; /* end marker */

  return data;
}

PRIVATE vrna_data_lin_t *
plist_to_ud_motif_prob(plist *pl, unsigned int length){

  int   i;
  plist *ptr;

  vrna_data_lin_t *data = NULL;

  data = (vrna_data_lin_t *)vrna_alloc(sizeof(vrna_data_lin_t) * (length + 1));

  for(ptr = pl; ptr->i > 0; ptr++){
    if(ptr->type == VRNA_PLIST_TYPE_UD_MOTIF){
      for(i = ptr->i; i <= ptr->j; i++){
        data[i - 1].value += ptr->p;
      }
    }
  }

  /*  go through the entire list, remove entries with 0 probability,
      and square-root probabilities again
  */
  unsigned int actual_length  = length;
  unsigned int actual_pos     = 1;
  for(i = 0; i < actual_length; i++, actual_pos++){
    if(data[i].value == 0.){
      memmove(data + i, data + i + 1, sizeof(vrna_data_lin_t) * (actual_length - i));
      actual_length--;
      i--;
    } else {
      data[i].position  = actual_pos;
      data[i].value     = sqrt(data[i].value);
      data[i].color.hue = 0.6;
      data[i].color.sat = 0.8;
      data[i].color.bri = 0.95;
    }
  }

  if(actual_length > 0){
    data[actual_length].position = 0; /* end marker */
    data = (vrna_data_lin_t *)vrna_realloc(data, sizeof(vrna_data_lin_t) * (actual_length + 1));
    return data;
  } else {
    free(data);
    return NULL;
  }
}


INLINE PRIVATE int
init_lin_data(lin_dat *dat)
{
  if (dat) {
    dat->north.size   = 0;
    dat->north.mem    = 8;
    dat->north.data   = (vrna_data_lin_t **)vrna_alloc(sizeof(vrna_data_lin_t *) * dat->north.mem);
    dat->north.titles = (char **)vrna_alloc(sizeof(char *) * dat->north.mem);

    dat->east.size   = 0;
    dat->east.mem    = 8;
    dat->east.data   = (vrna_data_lin_t **)vrna_alloc(sizeof(vrna_data_lin_t *) * dat->east.mem);
    dat->east.titles = (char **)vrna_alloc(sizeof(char *) * dat->east.mem);

    dat->south.size   = 0;
    dat->south.mem    = 8;
    dat->south.data   = (vrna_data_lin_t **)vrna_alloc(sizeof(vrna_data_lin_t *) * dat->south.mem);
    dat->south.titles = (char **)vrna_alloc(sizeof(char *) * dat->south.mem);

    dat->west.size   = 0;
    dat->west.mem    = 8;
    dat->west.data   = (vrna_data_lin_t **)vrna_alloc(sizeof(vrna_data_lin_t *) * dat->west.mem);
    dat->west.titles = (char **)vrna_alloc(sizeof(char *) * dat->west.mem);

    if ((dat->north.data == NULL) ||
        (dat->east.data == NULL) ||
        (dat->south.data == NULL) ||
        (dat->west.data == NULL) ||
        (dat->north.titles == NULL) ||
        (dat->east.titles == NULL) ||
        (dat->south.titles == NULL) ||
        (dat->west.titles == NULL))  {
      free(dat->north.data);
      free(dat->east.data);
      free(dat->south.data);
      free(dat->west.data);
      free(dat->north.titles);
      free(dat->east.titles);
      free(dat->south.titles);
      free(dat->west.titles);
      return 0; 
    }
  }

  return 1;
}


INLINE PRIVATE void
free_lin_data(lin_dat *dat)
{
  if (dat) {
      free(dat->north.data);
      free(dat->east.data);
      free(dat->south.data);
      free(dat->west.data);
      free(dat->north.titles);
      free(dat->east.titles);
      free(dat->south.titles);
      free(dat->west.titles);
  }
}


INLINE PRIVATE int
push_lin_data(lin_data_container  *c,
              vrna_data_lin_t     *data,
              char                *title)
{
  c->data[c->size]    = data;
  c->titles[c->size]  = title;

  if (++(c->size) == c->mem) {
    c->mem += 8;
    c->data = (vrna_data_lin_t **)vrna_realloc(c->data,
                                               sizeof(vrna_data_lin_t *) *
                                               c->mem);
    c->titles = (char **)vrna_realloc(c->titles,
                                      sizeof(char *) *
                                      c->mem);
  }

  /* check whether we successfully re-allocated memory */
  if ((c->data == NULL) ||
      (c->titles == NULL)) {
    free(c->data);
    free(c->titles);
    c->size = 0;
    c->mem  = 0;
    return 0;
  }

  return 1;
}


PRIVATE FILE *
PS_dot_common(const char    *seq,
              unsigned int  *nicks,
              const char    *wastlfile,
              char          *comment,
              int           winsize,
              unsigned int  options)
{

  /* write PS header etc for all dot plot variants */
  FILE *wastl;
  char *name, *c;

  wastl = fopen(wastlfile,"w");
  if (wastl==NULL) {
    vrna_log_warning("can't open %s for dot plot", wastlfile);
    return NULL; /* return 0 for failure */
  }
  name = strdup(wastlfile);
  if((c=strrchr(name, '_')))
    *c='\0';

  int bbox[4];
  if(winsize > 0){
    bbox[0] = 66;
    bbox[1] = 530;
    bbox[2] = 520;
    bbox[3] = 650;
  } else {
    bbox[0] = 66;
    bbox[1] = 211;
    bbox[2] = 518;
    bbox[3] = 662;
  }

  EPS_print_header(wastl, bbox, comment, options);

  EPS_print_title(wastl, name);

  print_PS_sequence(wastl, seq);

  if (winsize>0)
    fprintf(wastl,"/winSize %d def\n",winsize);

  if (nicks) {
    /* backward compatibility */
    fprintf(wastl,"/cutpoint %d def\n\n", nicks[0]);

    size_t cnt = 0;

    fprintf(wastl, "/nicks [ ");

    while (nicks[cnt])
      fprintf(wastl, "%d ", nicks[cnt++]);

    fprintf(wastl, "] def\n");
  }


  if (winsize>0)
  fprintf(wastl,"292 416 translate\n"
          "72 6 mul len 1 add winSize add 2 sqrt mul div dup scale\n");
  else
    fprintf(wastl,"72 216 translate\n"
          "72 6 mul len 1 add div dup scale\n");
  fprintf(wastl, "/Helvetica findfont 0.95 scalefont setfont\n\n");

  if (winsize>0) {
    fprintf(wastl, "%s", PS_dot_plot_macro_turn);
    fprintf(wastl,"0.5 dup translate\n"
          "drawseq_turn\n"
          "45 rotate\n\n");
  }
  else
    fprintf(wastl,"drawseq\n");

  free(name);
  return(wastl);
}

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*###########################################*/
/*# deprecated functions below              #*/
/*###########################################*/

int PS_dot_plot(char *string, char *wastlfile) {
  /* this is just a wrapper to call PS_dot_plot_list */
  int i, j, k, length, maxl, mf_num;
  plist *pl;
  plist *mf;

  i = 0;

  if ((string) &&
      (wastlfile) &&
      (pr) &&
      (iindx)) {
    length = strlen(string);
    maxl = 2*length;
    pl = (plist *)vrna_alloc(maxl*sizeof(plist));
    k=0;
    /*make plist out of pr array*/
    for (i=1; i<length; i++)
      for (j=i+1; j<=length; j++) {
        if (pr[iindx[i]-j]<PMIN) continue;
        if (k>=maxl-1) {
          maxl *= 2;
          pl = (plist *)vrna_realloc(pl,maxl*sizeof(plist));
        }
        pl[k].i = i;
        pl[k].j = j;
        pl[k].p = pr[iindx[i]-j];
        pl[k++].type = VRNA_PLIST_TYPE_BASEPAIR;
      }
    pl[k].i=0;
    pl[k].j=0;
    pl[k].p=0.;
    pl[k++].type=VRNA_PLIST_TYPE_BASEPAIR;
    /*make plist out of base_pair array*/
    mf_num = base_pair ? base_pair[0].i : 0;
    if (mf_num > 0) {
      mf = (plist *)vrna_alloc((mf_num+1)*sizeof(plist));
      for (k=0; k<mf_num; k++) {
        mf[k].i = base_pair[k+1].i;
        mf[k].j = base_pair[k+1].j;
        mf[k].p = 0.95*0.95;
        mf[k].type = VRNA_PLIST_TYPE_BASEPAIR;
      }
      mf[k].i=0;
      mf[k].j=0;
      mf[k].p=0.;
      mf[k].type=0;
    } else {
      mf = NULL;
    }
    i = vrna_plot_dp_PS_list(string, cut_point, wastlfile, pl, mf, "");
    free(mf);
    free(pl);
  }

  return (i);
}

#endif

