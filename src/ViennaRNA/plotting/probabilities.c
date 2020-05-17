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
#include "ViennaRNA/utils/alignments.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/plotting/probabilities.h"

#include "ViennaRNA/static/templates_postscript.h"

#include "ViennaRNA/plotting/ps_helpers.inc"

/*
#################################
# PRIVATE MACROS                #
#################################
*/

#define SIZE 452.
#define PMIN 0.00001

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
PRIVATE void  EPS_print_linear_data_top(FILE *eps, const char **data_title, vrna_data_lin_t **data);
PRIVATE void  EPS_print_linear_data_left(FILE *eps, const char **data_title, vrna_data_lin_t **data);
PRIVATE void  EPS_print_linear_data_bottom(FILE *eps, const char **data_title, vrna_data_lin_t **data);
PRIVATE void  EPS_print_linear_data_right(FILE *eps, const char **data_title, vrna_data_lin_t **data);
PRIVATE void  EPS_print_linear_data(FILE *eps, const char *varname, const char **data_title, vrna_data_lin_t **data);
PRIVATE vrna_data_lin_t *plist_to_accessibility(plist *pl, unsigned int length);
PRIVATE vrna_data_lin_t *plist_to_ud_motif_prob(plist *pl, unsigned int length);
PRIVATE void  EPS_print_sd_data(FILE *eps, vrna_ep_t *pl, vrna_ep_t *mf);

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
    } else if (pi[i].type == VRNA_PLIST_TYPE_BASEPAIR) {
      fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.6f ubox\n",
             pi[i].hue, pi[i].sat, pi[i].i, pi[i].j, sqrt(pi[i].p));

      if (pi[i].mfe)
        fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.4f lbox\n",
               pi[i].hue, pi[i].sat, pi[i].i, pi[i].j, pi[i].p);
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
  char          *seq_plain, *tmp, **seqs;
  unsigned int  *nicks, curr_length;
  int           pl_size, gq_num;
  plist         *pl1;

  nicks     = NULL;
  seq_plain = tmp = NULL;
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
      free(tmp);
      free(seqs[1]);

      /* add all remaining sequences (so far, we do not store the individual nicks for the remaining sequences */
      cnt = 2;
      while(seqs[cnt]) {
        nicks[cnt - 1] = curr_length + 1;
        vrna_strcat_printf(&seq_plain, "%s", seqs[cnt]);
        curr_length += strlen(seqs[cnt]);
        free(tmp);
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
  for(gq_num = pl_size = 0, pl1 = pl; pl1->i > 0; pl1++, pl_size++)
    if (pl1->type == VRNA_PLIST_TYPE_GQUAD) gq_num++;
  qsort(pl, pl_size, sizeof(plist), sort_plist_by_type_desc);
  /* sort all gquad triangles by probability to bring lower probs to the front */
  qsort(pl, gq_num, sizeof(plist), sort_plist_by_prob_asc);

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

  fh = fopen(filename, "w");
  if(!fh){
    vrna_message_warning("can't open %s for dot plot", filename);
    return 0; /* return 0 for failure */
  }

  comment = title = NULL;

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
      dp_add_lindata(lintop, ud_lin, lintoptitle, "Protein binding", lintop_num, lintop_size);
      dp_add_lindata(linleft, ud_lin, linlefttitle, "Protein binding", linleft_num, linleft_size);
      dp_add_lindata(linbottom, ud_lin, linbottomtitle, "Protein binding", linbottom_num, linbottom_size);
      dp_add_lindata(linright, ud_lin, linrighttitle, "Protein binding", linright_num, linright_size);
    }
  }

  if(options & VRNA_PLOT_PROBABILITIES_ACC){
    access = plist_to_accessibility(upper, strlen(sequence));
    dp_add_lindata(lintop, access, lintoptitle, "Accessibility", lintop_num, lintop_size);
  }

  if(auxdata){
    if(auxdata->top){
      for(i = 0; auxdata->top[i]; i++){
        dp_add_lindata(lintop, auxdata->top[i], lintoptitle, auxdata->top_title[i], lintop_num, lintop_size);
      }
    }
    if(auxdata->bottom){
      for(i = 0; auxdata->bottom[i]; i++){
        dp_add_lindata(linbottom, auxdata->bottom[i], linbottomtitle, auxdata->bottom_title[i], linbottom_num, linbottom_size);
      }
    }
    if(auxdata->left){
      for(i = 0; auxdata->left[i]; i++){
        dp_add_lindata(linleft, auxdata->left[i], linlefttitle, auxdata->left_title[i], linleft_num, linleft_size);
      }
    }
    if(auxdata->right){
      for(i = 0; auxdata->right[i]; i++){
        dp_add_lindata(linright, auxdata->right[i], linrighttitle, auxdata->right_title[i], linright_num, linright_size);
      }
    }
  }

  dp_finalize_lindata(lintop, lintoptitle, lintop_num);
  dp_finalize_lindata(linbottom, linbottomtitle, linbottom_num);
  dp_finalize_lindata(linleft, linlefttitle, linleft_num);
  dp_finalize_lindata(linright, linrighttitle, linright_num);

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
  EPS_print_linear_data_top(fh, (const char **)lintoptitle, lintop);
  EPS_print_linear_data_left(fh, (const char **)linlefttitle, linleft);
  EPS_print_linear_data_bottom(fh, (const char **)linbottomtitle, linbottom);
  EPS_print_linear_data_right(fh, (const char **)linrighttitle, linright);
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
             pi[i].hue, pi[i].sat, pi[i].i, pi[i].j, sqrt(pi[i].p));/*sqrt??*/

     if (pi[i].mfe)
       fprintf(wastl,"%1.2f %1.2f hsb %d %d %1.4f lbox\n",
               pi[i].hue, pi[i].sat, pi[i].i, pi[i].j, pi[i].p);
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
              pl[i].i, pl[i].j, sqrt(pl[i].p));
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
  for(gq_num = pl_size = 0, pl1 = pl; pl1->i > 0; pl1++, pl_size++)
    if(pl1->type == VRNA_PLIST_TYPE_GQUAD) gq_num++;

  qsort(pl, pl_size, sizeof(plist), sort_plist_by_type_desc);

  /* sort all gquad triangles by probability to bring lower probs to the front */
  qsort(pl, gq_num, sizeof(plist), sort_plist_by_prob_asc);

  /* print triangles for g-quadruplexes in upper half */
  fprintf(eps,"\n%%start of quadruplex data\n");

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
  for (pl1=pl; pl1->i > 0; pl1++) {
    if(pl1->type == VRNA_PLIST_TYPE_H_MOTIF){
      tmp = sqrt(pl1->p);
      fprintf(eps, "%d %d %1.9f uHmotif\n", pl1->i, pl1->j, tmp);
    }
  }
  for (pl1=mf; pl1->i > 0; pl1++) {
    if(pl1->type == VRNA_PLIST_TYPE_H_MOTIF){
      tmp = sqrt(pl1->p);
      fprintf(eps, "%d %d %1.9f lHmotif\n", pl1->i, pl1->j, tmp);
    }
  }

  /* print triangles for interior loop motifs in upper half */
  fprintf(eps,"\n%%start of Imotif data\n");
  int   a,b;
  float ppp;
  a = b = 0;
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
  for (pl1 = pl; pl1->i>0; pl1++) {
    tmp = sqrt(pl1->p);
    if(pl1->type == VRNA_PLIST_TYPE_BASEPAIR)
        fprintf(eps,"%d %d %1.9f ubox\n", pl1->i, pl1->j, tmp);
  }


  /* print boxes in lower left half (mfe) */
  for (pl1=mf; pl1->i>0; pl1++) {
    tmp = sqrt(pl1->p);
    if(pl1->type == VRNA_PLIST_TYPE_BASEPAIR)
      fprintf(eps,"%d %d %1.7f lbox\n", pl1->i, pl1->j, tmp);
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
  for(pl1 = pl; pl1->i > 0; pl1++){
    if(pl1->type == VRNA_PLIST_TYPE_UD_MOTIF){
      tmp = sqrt(pl1->p);
      fprintf(eps, "%d %d %1.9f uUDmotif\n", pl1->i, pl1->j, tmp);
    }
  }
  for(pl1 = mf; pl1->i > 0; pl1++){
    if(pl1->type == VRNA_PLIST_TYPE_UD_MOTIF){
      tmp = sqrt(pl1->p);
      fprintf(eps, "%d %d %1.9f lUDmotif\n", pl1->i, pl1->j, tmp);
    }
  }
}


PRIVATE void
EPS_print_linear_data_top(FILE            *eps,
                          const char      **data_title,
                          vrna_data_lin_t **data){

  EPS_print_linear_data(eps, "topData", data_title, data);
}


PRIVATE void
EPS_print_linear_data_left( FILE            *eps,
                            const char      **data_title,
                            vrna_data_lin_t **data){

  EPS_print_linear_data(eps, "leftData", data_title, data);
}


PRIVATE void
EPS_print_linear_data_bottom( FILE            *eps,
                              const char      **data_title,
                              vrna_data_lin_t **data){

  EPS_print_linear_data(eps, "bottomData", data_title, data);
}


PRIVATE void
EPS_print_linear_data_right(FILE            *eps,
                            const char      **data_title,
                            vrna_data_lin_t **data){

  EPS_print_linear_data(eps, "rightData", data_title, data);
}


PRIVATE void
EPS_print_linear_data(FILE            *eps,
                      const char      *varname,
                      const char      **data_title,
                      vrna_data_lin_t **data){

  int             n, i;
  vrna_data_lin_t *ptr;

  /* count number of data sets */
  for(n = 0; data_title[n]; n++);

  fprintf(eps, "/%s [\n", varname);
  for(i = 0; i < n; i++){
    fprintf(eps, "[ (%s)\n", data_title[i]);
    for(ptr = data[i]; ptr->position > 0; ptr++){
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
    vrna_message_warning("can't open %s for dot plot", wastlfile);
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
      pl[k++].p = pr[iindx[i]-j];
    }
  pl[k].i=0;
  pl[k].j=0;
  pl[k++].p=0.;
  /*make plist out of base_pair array*/
  mf_num = base_pair ? base_pair[0].i : 0;
  mf = (plist *)vrna_alloc((mf_num+1)*sizeof(plist));
  for (k=0; k<mf_num; k++) {
    mf[k].i = base_pair[k+1].i;
    mf[k].j = base_pair[k+1].j;
    mf[k].p = 0.95*0.95;
  }
  mf[k].i=0;
  mf[k].j=0;
  mf[k].p=0.;
  i = PS_dot_plot_list(string, wastlfile, pl, mf, "");
  free(mf);
  free(pl);
  return (i);
}

#endif

