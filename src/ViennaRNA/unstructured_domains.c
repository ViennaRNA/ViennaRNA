/** \file ligands_up.c **/

/*
                  Ligand binding to unpaired stretches

                  This file contains everything necessary to
                  deal with the default implementation for ligand
                  binding to unpaired stretches of an RNA secondary
                  structure

                  c 2016 Ronny Lorenz

                  ViennaRNA package
*/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "ViennaRNA/utils.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/unstructured_domains.h"

/*
#################################
# PRIVATE MACROS                #
#################################
*/

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/

/*
#################################
# PRIVATE VARIABLES/STRUCTS     #
#################################
*/

/*
 *  Default data structure for ligand binding to unpaired stretches
 */
struct ligands_up_data_default {

  /*
  **********************************
    pre-computed position-wise
    motif list
  **********************************
  */
  int           n;
  int           **motif_list_ext;
  int           **motif_list_hp;
  int           **motif_list_int;
  int           **motif_list_mb;

  int           *dG;
  int           *len;

  /*
  **********************************
    below are DP matrices to store
    the production rule results
  **********************************
  */
  int           *energies_ext;
  int           *energies_hp;
  int           *energies_int;
  int           *energies_mb;
  FLT_OR_DBL    *exp_energies_ext;
  FLT_OR_DBL    *exp_energies_hp;
  FLT_OR_DBL    *exp_energies_int;
  FLT_OR_DBL    *exp_energies_mb;

};

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/

PRIVATE void        remove_ligands_up(vrna_fold_compound_t *vc);
PRIVATE void        init_ligands_up(vrna_fold_compound_t *vc);

PRIVATE void        add_ligand_motif(vrna_fold_compound_t *vc, const char *motif, double motif_en, unsigned int loop_type);
PRIVATE void        remove_default_data(void *d);

PRIVATE void        default_prod_rule(vrna_fold_compound_t *vc, void *d);
PRIVATE int         default_energy(vrna_fold_compound_t *vc, int i, int j, unsigned int loop_type, void *d);
PRIVATE int         default_energy_ext_motif(int i, int j, struct ligands_up_data_default *data);
PRIVATE int         default_energy_hp_motif(int i, int j, struct ligands_up_data_default *data);
PRIVATE int         default_energy_int_motif(int i, int j, struct ligands_up_data_default *data);
PRIVATE int         default_energy_mb_motif(int i, int j, struct ligands_up_data_default *data);

PRIVATE void        default_exp_prod_rule(vrna_fold_compound_t *vc, void *d);
PRIVATE FLT_OR_DBL  default_exp_energy( vrna_fold_compound_t *vc, int i, int j, unsigned int loop_type, void *d);
PRIVATE FLT_OR_DBL  default_exp_energy_ext_motif(int i, int j, struct ligands_up_data_default *data);
PRIVATE FLT_OR_DBL  default_exp_energy_hp_motif(int i, int j, struct ligands_up_data_default *data);
PRIVATE FLT_OR_DBL  default_exp_energy_int_motif(int i, int j, struct ligands_up_data_default *data);
PRIVATE FLT_OR_DBL  default_exp_energy_mb_motif(int i, int j, struct ligands_up_data_default *data);

PRIVATE int         *get_motifs(vrna_fold_compound_t *vc, int i, unsigned int loop_type);
PRIVATE void        free_default_data_matrices(struct ligands_up_data_default *data);
PRIVATE void        prepare_matrices( vrna_fold_compound_t *vc, struct ligands_up_data_default *data);
PRIVATE struct ligands_up_data_default *get_default_data(void);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/
PUBLIC void
vrna_ud_remove( vrna_fold_compound_t *vc){

  if(vc && vc->domains_up)
    remove_ligands_up(vc);
}

PUBLIC void
vrna_ud_set_data( vrna_fold_compound_t  *vc,
                          void *data,
                          vrna_callback_free_auxdata  *free){

  if(vc){
    /* init if not already present */
    if(!vc->domains_up)
      init_ligands_up(vc);

    /* free previous data if 'free_data' function present */
    if(vc->domains_up->free_data)
      vc->domains_up->free_data(vc->domains_up->data);

    /* set new data and free callback */
    vc->domains_up->free_data = free;
    vc->domains_up->data      = data;
  }
}

PUBLIC void
vrna_ud_set_prod_rule(vrna_fold_compound_t  *vc,
                              vrna_callback_ud_production *rule){

  if(vc){
    /* init if not already present */
    if(!vc->domains_up)
      init_ligands_up(vc);

    /* set new callback */
    vc->domains_up->prod_cb = rule;
  }
}

PUBLIC void
vrna_ud_set_exp_prod_rule(vrna_fold_compound_t  *vc,
                                  vrna_callback_ud_exp_production *rule){

  if(vc){
    /* init if not already present */
    if(!vc->domains_up)
      init_ligands_up(vc);

    /* set new callback */
    vc->domains_up->exp_prod_cb = rule;
  }
}

PUBLIC void
vrna_ud_set_energy( vrna_fold_compound_t *vc,
                            vrna_callback_ud_energy *e){

  if(vc){
    /* init if not already present */
    if(!vc->domains_up)
      init_ligands_up(vc);

    /* set new callback */
    vc->domains_up->energy_cb = e;
  }
}

PUBLIC void
vrna_ud_set_exp_energy( vrna_fold_compound_t *vc,
                                vrna_callback_ud_exp_energy *exp_e){

  if(vc){
    /* init if not already present */
    if(!vc->domains_up)
      init_ligands_up(vc);

    /* set new callback */
    vc->domains_up->exp_energy_cb = exp_e;
  }
}

PUBLIC void
vrna_ud_add_motif(vrna_fold_compound_t*vc,
                  const char *motif,
                  double motif_en,
                  unsigned int loop_type){

  if(vc){
    if(!vc->domains_up){
      init_ligands_up(vc);
      /* set all callbacks and stuff to default settings */
      vc->domains_up->prod_cb       = &default_prod_rule;
      vc->domains_up->exp_prod_cb   = &default_exp_prod_rule;
      vc->domains_up->energy_cb     = &default_energy;
      vc->domains_up->exp_energy_cb = &default_exp_energy;
      vc->domains_up->data          = get_default_data();
      vc->domains_up->free_data     = &remove_default_data;
    }
    add_ligand_motif(vc, motif, motif_en, loop_type);
  }
}
/*
#####################################
# BEGIN OF STATIC HELPER FUNCTIONS  #
#####################################
*/
PRIVATE struct ligands_up_data_default *
get_default_data(void){

  struct ligands_up_data_default *data = vrna_alloc(sizeof(struct ligands_up_data_default));
  data->n                 = 0;
  data->motif_list_ext    = NULL;
  data->motif_list_hp     = NULL;
  data->motif_list_int    = NULL;
  data->motif_list_mb     = NULL;
  data->dG                = NULL;
  data->energies_ext      = NULL;
  data->energies_hp       = NULL;
  data->energies_int      = NULL;
  data->energies_mb       = NULL;
  data->exp_energies_ext  = NULL;
  data->exp_energies_hp   = NULL;
  data->exp_energies_int  = NULL;
  data->exp_energies_mb   = NULL;

  return data;
}

PRIVATE void
remove_ligands_up(vrna_fold_compound_t *vc){

  int i;

  /* free auxiliary data */
  if(vc->domains_up->free_data)
    vc->domains_up->free_data(vc->domains_up->data);

  for( i = 0; i < vc->domains_up->motif_count; i++ ){
    free(vc->domains_up->motif[i]);
  }
  free(vc->domains_up->motif);
  free(vc->domains_up->motif_size);
  free(vc->domains_up->motif_en);
  free(vc->domains_up->motif_type);

  free(vc->domains_up);

  vc->domains_up = NULL;
}

PRIVATE void
init_ligands_up(vrna_fold_compound_t *vc){

  vc->domains_up = (vrna_ud_t *)vrna_alloc(sizeof(vrna_ud_t));

  vc->domains_up->motif_count   = 0;
  vc->domains_up->motif         = NULL;
  vc->domains_up->motif_size    = NULL;
  vc->domains_up->motif_en      = NULL;
  vc->domains_up->motif_type    = NULL;
  vc->domains_up->prod_cb       = NULL;
  vc->domains_up->exp_prod_cb   = NULL;
  vc->domains_up->energy_cb     = NULL;
  vc->domains_up->exp_energy_cb = NULL;
  vc->domains_up->data          = NULL;
  vc->domains_up->free_data     = NULL;
}

/*
**********************************
  Default implementation for
  ligand binding to unpaired
  stretches follows below
**********************************
*/

PRIVATE void
add_ligand_motif( vrna_fold_compound_t *vc,
                  const char *motif,
                  double motif_en,
                  unsigned int loop_type){

  vc->domains_up->motif                                   = (char **)vrna_realloc(vc->domains_up->motif, sizeof(char *) * (vc->domains_up->motif_count + 1));
  vc->domains_up->motif[vc->domains_up->motif_count]      = strdup(motif);
  vc->domains_up->motif_size                              = (unsigned int *)vrna_realloc(vc->domains_up->motif_size, sizeof(unsigned int *) * (vc->domains_up->motif_count + 1));
  vc->domains_up->motif_size[vc->domains_up->motif_count] = strlen(motif);
  vc->domains_up->motif_en                                = (double *)vrna_realloc(vc->domains_up->motif_en, sizeof(double) * (vc->domains_up->motif_count + 1));
  vc->domains_up->motif_en[vc->domains_up->motif_count]   = motif_en;
  vc->domains_up->motif_type                              = (unsigned int *)vrna_realloc(vc->domains_up->motif_type, sizeof(double) * (vc->domains_up->motif_count + 1));
  vc->domains_up->motif_type[vc->domains_up->motif_count] = loop_type;
  vc->domains_up->motif_count++;
}

PRIVATE void
remove_default_data(void *d){

  int i;
  struct ligands_up_data_default  *data;

  data = (struct ligands_up_data_default *)d;

  if(data->motif_list_ext){
    for(i=0; i <= data->n; i++)
      free(data->motif_list_ext[i]);
    free(data->motif_list_ext);
  }
  if(data->motif_list_hp){
    for(i=0; i <= data->n; i++)
      free(data->motif_list_hp[i]);
    free(data->motif_list_hp);
  }
  if(data->motif_list_int){
    for(i=0; i <= data->n; i++)
      free(data->motif_list_int[i]);
    free(data->motif_list_int);
  }
  if(data->motif_list_mb){
    for(i=0; i <= data->n; i++)
      free(data->motif_list_mb[i]);
    free(data->motif_list_mb);
  }

  free_default_data_matrices(data);

  free(data->dG);
  free(data->len);
}

PRIVATE void
free_default_data_matrices(struct ligands_up_data_default *data){

  /* the following four pointers may point to the same memory */
  if(data->energies_ext){
    /* check whether one of the other b* points to the same memory location */
    if(data->energies_ext == data->energies_hp)
      data->energies_hp = NULL;
    if(data->energies_ext == data->energies_int)
      data->energies_int = NULL;
    if(data->energies_ext == data->energies_mb)
      data->energies_mb = NULL;
    free(data->energies_ext);
    data->energies_ext = NULL;
  }
  if(data->energies_hp){
    /* check whether one of the other b* points to the same memory location */
    if(data->energies_hp == data->energies_int)
      data->energies_int = NULL;
    if(data->energies_hp == data->energies_mb)
      data->energies_mb = NULL;
    free(data->energies_hp);
    data->energies_hp = NULL;
  }
  if(data->energies_int){
    /* check whether one of the other b* points to the same memory location */
    if(data->energies_int == data->energies_mb)
      data->energies_mb = NULL;
    free(data->energies_int);
    data->energies_int = NULL;
  }
  free(data->energies_mb);
  data->energies_mb = NULL;

  /* the following four pointers may point to the same memory */
  if(data->exp_energies_ext){
    /* check whether one of the other b* points to the same memory location */
    if(data->exp_energies_ext == data->exp_energies_hp)
      data->exp_energies_hp = NULL;
    if(data->exp_energies_ext == data->exp_energies_int)
      data->exp_energies_int = NULL;
    if(data->exp_energies_ext == data->exp_energies_mb)
      data->exp_energies_mb = NULL;
    free(data->exp_energies_ext);
  }
  if(data->exp_energies_hp){
    /* check whether one of the other b* points to the same memory location */
    if(data->exp_energies_hp == data->exp_energies_int)
      data->exp_energies_int = NULL;
    if(data->exp_energies_hp == data->exp_energies_mb)
      data->exp_energies_mb = NULL;
    free(data->exp_energies_hp);
  }
  if(data->exp_energies_int){
    /* check whether one of the other b* points to the same memory location */
    if(data->exp_energies_int == data->exp_energies_mb)
      data->exp_energies_mb = NULL;
    free(data->exp_energies_int);
  }
  free(data->exp_energies_mb);
}

PRIVATE int *
get_motifs(vrna_fold_compound_t *vc, int i, unsigned int loop_type){

  int               k, j, u, n, *motif_list, cnt, guess;
  char              *sequence;
  vrna_ud_t *ligands_up;

  sequence    = vc->sequence;
  n           = (int)vc->length;
  ligands_up  = vc->domains_up;

  cnt         = 0;
  guess       = ligands_up->motif_count;
  motif_list  = (int *)vrna_alloc(sizeof(int) * (guess + 1));

  /* collect list of motif numbers we find that start at position i */
  for(k = 0; k < ligands_up->motif_count; k++){

    if(!(ligands_up->motif_type[k] & loop_type))
      continue;

    j = i + ligands_up->motif_size[k] - 1;
    if(j <= n){ /* only consider motif that does not exceed sequence length (does not work for circular RNAs!) */
      for(u = i; u <= j; u++){
        if(!vrna_nucleotide_IUPAC_identity(sequence[u-1], ligands_up->motif[k][u-i]))
          break;
      }
      if(u > j) /* got a complete motif match */
        motif_list[cnt++] = k;
    }
  }

  if(cnt == 0){
    free(motif_list);
    return NULL;
  }
  
  motif_list = (int *)vrna_realloc(motif_list, sizeof(int) * (cnt + 1));
  motif_list[cnt] = -1; /* end of list marker */

  return motif_list;
}

PRIVATE void
prepare_matrices( vrna_fold_compound_t *vc,
                  struct ligands_up_data_default *data){

  int                             i,j,k,n,size;
  vrna_ud_t               *ligands_up;

  n             = (int)vc->length;
  size          = ((n+1)*(n+2))/2 + 1;
  ligands_up    = vc->domains_up;

  free_default_data_matrices(data);

  /* here we save memory by re-using DP matrices */
  unsigned int lt[4] = {  VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                          VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                          VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                          VRNA_UNSTRUCTURED_DOMAIN_ML_LOOP };
  int **m[4], *mx;
  m[0] = &data->energies_ext;
  m[1] = &data->energies_hp;
  m[2] = &data->energies_int;
  m[3] = &data->energies_mb;

  for(i=0; i<4; i++){
    unsigned int *col,*col2;
    if(*(m[i]))
      continue;

    mx      = (int *)vrna_alloc(sizeof(int) * size);
    col     = (unsigned int *)vrna_alloc(sizeof(unsigned int) * ligands_up->motif_count);
    col2    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * ligands_up->motif_count);
    *(m[i]) = mx;

    for(k = 0; k < ligands_up->motif_count; k++)
      col[k] = ligands_up->motif_type[k] & lt[i];

    /* check if any of the remaining DP matrices can point to the same location */
    for(j=i+1;j<4;j++){
      for(k = 0; k < ligands_up->motif_count; k++){
        col2[k] = ligands_up->motif_type[k] & lt[j];
        if(col[k] != col2[k])
          break;
      }
      if(k == ligands_up->motif_count){
        *(m[j]) = mx;
      }
    }

    free(col);
    free(col2);
  }
}

PRIVATE void
default_prod_rule(vrna_fold_compound_t *vc,
                  void *d){

  int                             i,j,k,l,u,n,size,e_ext, e_hp, e_int, e_mb,en,en2,*idx;
  unsigned int                    loop_type;
  vrna_ud_t               *ligands_up;
  struct ligands_up_data_default  *data;

  int           *energies_ext;
  int           *energies_hp;
  int           *energies_int;
  int           *energies_mb;

  n             = (int)vc->length;
  size          = ((n+1)*(n+2))/2 + 1;
  idx           = vc->jindx;
  ligands_up    = vc->domains_up;
  data          = (struct ligands_up_data_default *)d;

  prepare_matrices(vc, data);

  energies_ext  = data->energies_ext;
  energies_hp   = data->energies_hp;
  energies_int  = data->energies_int;
  energies_mb   = data->energies_mb;

  data->n = n;
  /*  create motif_list for associating a nucleotide position with all
      motifs that start there
  */
  data->motif_list_ext  = (int **)vrna_alloc(sizeof(int *) * (n+1));
  data->motif_list_hp   = (int **)vrna_alloc(sizeof(int *) * (n+1));
  data->motif_list_int  = (int **)vrna_alloc(sizeof(int *) * (n+1));
  data->motif_list_mb   = (int **)vrna_alloc(sizeof(int *) * (n+1));
  data->motif_list_ext[0] = NULL;
  data->motif_list_hp[0]  = NULL;
  data->motif_list_int[0] = NULL;
  data->motif_list_mb[0]  = NULL;
  for(i = 1; i <= n; i++){
    data->motif_list_ext[i] = get_motifs(vc, i, VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP);
    data->motif_list_hp[i]  = get_motifs(vc, i, VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP);
    data->motif_list_int[i] = get_motifs(vc, i, VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP);
    data->motif_list_mb[i]  = get_motifs(vc, i, VRNA_UNSTRUCTURED_DOMAIN_ML_LOOP);
  }

  /*  precompute energy contributions of the motifs */
  data->dG = (int *)vrna_alloc(sizeof(int) * ligands_up->motif_count);
  for(i = 0; i < ligands_up->motif_count; i++)
    data->dG[i] = (int)roundf(ligands_up->motif_en[i] * 100.);

  /*  store length of motifs in 'data' */
  data->len = (int *)vrna_alloc(sizeof(int) * ligands_up->motif_count);
  for(i = 0; i < ligands_up->motif_count; i++)
    data->len[i] = ligands_up->motif_size[i];

  /* now we can start to fill the DP matrices */
  for(i=n; i>0; i--){
    int *list_ext = data->motif_list_ext[i];
    int *list_hp  = data->motif_list_hp[i];
    int *list_int = data->motif_list_int[i];
    int *list_mb  = data->motif_list_mb[i];
    for(j=i;j<=n;j++){
      if(i < n){
        e_ext = energies_ext[idx[j]+i+1];
        e_hp  = energies_hp[idx[j]+i+1];
        e_int = energies_int[idx[j]+i+1];
        e_mb  = energies_mb[idx[j]+i+1];
      } else {
        e_ext = 0;
        e_hp  = 0;
        e_int = 0;
        e_mb  = 0;
      }
      if(list_ext){
        for(k = 0; -1 != (l = list_ext[k]); k++){
          u   = i + data->len[l] - 1;
          en  = data->dG[l];
          if(u < j){
            en2 = en + energies_ext[idx[j]+u+1];
            e_ext = MIN2(e_ext, en2);
          } else if(u == j){
            e_ext = MIN2(e_ext, en);
          }
        }
      }
      if(list_hp){
        for(k = 0; -1 != (l = list_hp[k]); k++){
          u   = i + data->len[l] - 1;
          en  = data->dG[l];
          if(u < j){
            en2 = en + energies_hp[idx[j]+u+1];
            e_hp = MIN2(e_hp, en2);
          } else if(u == j){
            e_hp = MIN2(e_hp, en);
          }
        }
      }
      if(list_int){
        for(k = 0; -1 != (l = list_int[k]); k++){
          u   = i + data->len[l] - 1;
          en  = data->dG[l];
          if(u < j){
            en2 = en + energies_int[idx[j]+u+1];
            e_int = MIN2(e_int, en2);
          } else if(u == j){
            e_int = MIN2(e_int, en);
          }
        }
      }
      if(list_mb){
        for(k = 0; -1 != (l = list_mb[k]); k++){
          u   = i + data->len[l] - 1;
          en  = data->dG[l];
          if(u < j){
            en2 = en + energies_mb[idx[j]+u+1];
            e_mb = MIN2(e_mb, en2);
          } else if(u == j){
            e_mb = MIN2(e_mb, en);
          }
        }
      }

      energies_ext[idx[j]+i]  = e_ext;
      energies_hp[idx[j]+i]   = e_hp;
      energies_int[idx[j]+i]  = e_int;
      energies_mb[idx[j]+i]   = e_mb;
    }
  }
}

PRIVATE void
default_exp_prod_rule(vrna_fold_compound_t *vc,
                      void *d){

  struct ligands_up_data_default *data = (struct ligands_up_data_default *)d;

}

PRIVATE int
default_energy( vrna_fold_compound_t *vc,
                int i,
                int j,
                unsigned int loop_type,
                void *d){

  int en, ij, *idx = vc->jindx;
  struct ligands_up_data_default *data = (struct ligands_up_data_default *)d;

  en = INF;
  ij = idx[j] + i;

  if(j < i)
    return 0;

  if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_MOTIF){
    switch(loop_type){
      case VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP:   en = default_energy_ext_motif(i, j, data);
                                                break;

      case VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP:    en = default_energy_hp_motif(i, j, data);
                                                break;

      case VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP:   en = default_energy_int_motif(i, j, data);
                                                break;

      case VRNA_UNSTRUCTURED_DOMAIN_ML_LOOP:    en = default_energy_mb_motif(i, j, data);
                                                break;
    }
  } else {
    switch(loop_type){
      case VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP:   if(data->energies_ext)
                                                  en = data->energies_ext[ij];
                                                break;

      case VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP:    if(data->energies_hp)
                                                  en = data->energies_hp[ij];
                                                break;

      case VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP:   if(data->energies_int)
                                                  en = data->energies_int[ij];
                                                break;

      case VRNA_UNSTRUCTURED_DOMAIN_ML_LOOP:    if(data->energies_mb)
                                                  en = data->energies_mb[ij];
                                                break;
    }
  }

  return en;
}

PRIVATE FLT_OR_DBL
default_exp_energy( vrna_fold_compound_t *vc,
                    int i,
                    int j,
                    unsigned int loop_type,
                    void *d){

  FLT_OR_DBL                      q;
  int                             ij, *idx;
  struct ligands_up_data_default  *data;

  q     = 0;
  data  = (struct ligands_up_data_default *)d;
  idx   = vc->iindx;
  ij    = idx[i] - j;

  if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_MOTIF){
    switch(loop_type){
      case VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP: q = default_exp_energy_ext_motif(i, j, data);
                                              break;

      case VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP:  q = default_exp_energy_hp_motif(i, j, data);
                                              break;

      case VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP: q = default_exp_energy_int_motif(i, j, data);
                                              break;

      case VRNA_UNSTRUCTURED_DOMAIN_ML_LOOP:  q = default_exp_energy_mb_motif(i, j, data);
                                              break;
    }
  } else {
    switch(loop_type){
      case VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP: if(data->exp_energies_ext)
                                              q = data->exp_energies_ext[ij];
                                              break;

      case VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP:  if(data->exp_energies_hp)
                                              q = data->exp_energies_hp[ij];
                                              break;

      case VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP: if(data->exp_energies_int)
                                              q = data->exp_energies_int[ij];
                                              break;

      case VRNA_UNSTRUCTURED_DOMAIN_ML_LOOP:  if(data->exp_energies_mb)
                                              q = data->exp_energies_mb[ij];
                                              break;
    }
  }

  return q;
}

PRIVATE int
default_energy_ext_motif( int i,
                          int j,
                          struct ligands_up_data_default *data){

  int k, m;
  int e = INF;

  if(data->motif_list_ext[i]){
    k = 0;
    while(-1 != (m = data->motif_list_ext[i][k])){
      if((i + data->len[m] - 1) == j)
        return data->dG[m];
      k++;
    }
  }
  return e;
}

PRIVATE int
default_energy_hp_motif(int i,
                        int j,
                        struct ligands_up_data_default *data){

  int k, m;
  int e = INF;

  if(data->motif_list_hp[i]){
    k = 0;
    while(-1 != (m = data->motif_list_hp[i][k])){
      if((i + data->len[m] - 1) == j)
        return data->dG[m];
      k++;
    }
  }
  return e;
}

PRIVATE int
default_energy_int_motif( int i,
                          int j,
                          struct ligands_up_data_default *data){

  int k, m;
  int e = INF;

  if(data->motif_list_int[i]){
    k = 0;
    while(-1 != (m = data->motif_list_int[i][k])){
      if((i + data->len[m] - 1) == j)
        return data->dG[m];
      k++;
    }
  }
  return e;
}

PRIVATE int
default_energy_mb_motif(int i,
                        int j,
                        struct ligands_up_data_default *data){

  int k, m;
  int e = INF;

  if(data->motif_list_mb[i]){
    k = 0;
    while(-1 != (m = data->motif_list_mb[i][k])){
      if((i + data->len[m] - 1) == j)
        return data->dG[m];
      k++;
    }
  }
  return e;
}

PRIVATE FLT_OR_DBL
default_exp_energy_ext_motif( int i,
                              int j,
                              struct ligands_up_data_default *data){

  FLT_OR_DBL q = 0;

  return q;
}

PRIVATE FLT_OR_DBL
default_exp_energy_hp_motif(int i,
                            int j,
                            struct ligands_up_data_default *data){

  FLT_OR_DBL q = 0;

  return q;
}

PRIVATE FLT_OR_DBL
default_exp_energy_int_motif( int i,
                              int j,
                              struct ligands_up_data_default *data){

  FLT_OR_DBL q = 0;

  return q;
}

PRIVATE FLT_OR_DBL
default_exp_energy_mb_motif(int i,
                            int j,
                            struct ligands_up_data_default *data){

  FLT_OR_DBL q = 0;

  return q;
}
