/** \file unstructured_domains.c **/

/*
                  Unstructured domains

                  This file contains everything necessary to
                  deal with the default implementation for unstructured
                  domains in secondary structures. This feature enables,
                  for instance, ligand binding to unpaired stretches of
                  an RNA secondary structure.

                  c 2016 Ronny Lorenz

                  ViennaRNA package
*/

#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
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

struct default_outside {
  int         motif_num;
  FLT_OR_DBL  exp_energy;
};

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
  FLT_OR_DBL    *exp_dG;
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

  /*
  **********************************
    below are lists to store the
    outside partition function for
    each motif starting at each
    position
  **********************************
  */
  unsigned int            *outside_ext_count;
  struct default_outside  **outside_ext;
  unsigned int            *outside_hp_count;
  struct default_outside  **outside_hp;
  unsigned int            *outside_int_count;
  struct default_outside  **outside_int;
  unsigned int            *outside_mb_count;
  struct default_outside  **outside_mb;

  FLT_OR_DBL  (*default_cb[32])(int, int, struct ligands_up_data_default *);
  FLT_OR_DBL  *exp_e_mx[32];
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

/* default implementations for unstructured domains feature */
PRIVATE void        default_prod_rule(vrna_fold_compound_t *vc, void *d);
PRIVATE void        default_exp_prod_rule(vrna_fold_compound_t *vc, void *d);
PRIVATE int         default_energy(vrna_fold_compound_t *vc, int i, int j, unsigned int loop_type, void *d);
PRIVATE FLT_OR_DBL  default_exp_energy( vrna_fold_compound_t *vc, int i, int j, unsigned int loop_type, void *d);
PRIVATE void        default_probs_add(vrna_fold_compound_t *vc, int i, int j, unsigned int loop_type, FLT_OR_DBL exp_energy, void *data);
PRIVATE FLT_OR_DBL  default_probs_get(vrna_fold_compound_t *vc, int i, int j, unsigned int loop_type, int motif, void *data);

/* helper functions for default implementatations of unstructured domains feature */
PRIVATE int         default_energy_ext_motif(int i, int j, struct ligands_up_data_default *data);
PRIVATE int         default_energy_hp_motif(int i, int j, struct ligands_up_data_default *data);
PRIVATE int         default_energy_int_motif(int i, int j, struct ligands_up_data_default *data);
PRIVATE int         default_energy_mb_motif(int i, int j, struct ligands_up_data_default *data);
PRIVATE FLT_OR_DBL  default_exp_energy_ext_motif(int i, int j, struct ligands_up_data_default *data);
PRIVATE FLT_OR_DBL  default_exp_energy_hp_motif(int i, int j, struct ligands_up_data_default *data);
PRIVATE FLT_OR_DBL  default_exp_energy_int_motif(int i, int j, struct ligands_up_data_default *data);
PRIVATE FLT_OR_DBL  default_exp_energy_mb_motif(int i, int j, struct ligands_up_data_default *data);
PRIVATE void        free_default_data_matrices(struct ligands_up_data_default *data);
PRIVATE void        free_default_data_exp_matrices(struct ligands_up_data_default *data);
PRIVATE void        prepare_matrices( vrna_fold_compound_t *vc, struct ligands_up_data_default *data);
PRIVATE void        prepare_exp_matrices( vrna_fold_compound_t *vc, struct ligands_up_data_default *data);
PRIVATE struct ligands_up_data_default *get_default_data(void);
PRIVATE void        prepare_default_data(vrna_fold_compound_t *vc, struct ligands_up_data_default *data);
PRIVATE void        free_default_data(struct ligands_up_data_default *data);
PRIVATE int         *get_motifs(vrna_fold_compound_t *vc, int i, unsigned int loop_type);
PRIVATE void        annotate_ud(vrna_fold_compound_t *vc, int start, int end, char l, vrna_ud_motif_t **list, int *list_size, int *list_pos);


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
vrna_ud_set_data( vrna_fold_compound_t                *vc,
                          void                        *data,
                          vrna_callback_free_auxdata  *free_cb){

  if(vc){
    /* init if not already present */
    if(!vc->domains_up)
      init_ligands_up(vc);

    /* free previous data if 'free_data' function present */
    if(vc->domains_up->free_data)
      vc->domains_up->free_data(vc->domains_up->data);

    /* set new data and free callback */
    vc->domains_up->free_data = free_cb;
    vc->domains_up->data      = data;
  }
}

PUBLIC void
vrna_ud_set_prod_rule_cb( vrna_fold_compound_t        *vc,
                          vrna_callback_ud_production *pre_cb,
                          vrna_callback_ud_energy     *e_cb){

  if(vc){
    /* init if not already present */
    if(!vc->domains_up)
      init_ligands_up(vc);

    /* set new callback */
    vc->domains_up->prod_cb   = pre_cb;
    vc->domains_up->energy_cb = e_cb;
  }
}

PUBLIC void
vrna_ud_set_exp_prod_rule_cb( vrna_fold_compound_t            *vc,
                              vrna_callback_ud_exp_production *pre_cb,
                              vrna_callback_ud_exp_energy     *exp_e_cb){

  if(vc){
    /* init if not already present */
    if(!vc->domains_up)
      init_ligands_up(vc);

    /* set new callback */
    vc->domains_up->exp_prod_cb   = pre_cb;
    vc->domains_up->exp_energy_cb = exp_e_cb;
  }
}


PUBLIC void
vrna_ud_set_prob_cb(vrna_fold_compound_t        *vc,
                    vrna_callback_ud_probs_add  *setter,
                    vrna_callback_ud_probs_get  *getter){

  if(vc){
    /* init if not already present */
    if(!vc->domains_up)
      init_ligands_up(vc);

    /* set new callback */
    vc->domains_up->probs_add = setter;
    vc->domains_up->probs_get = getter;
  }
}

PUBLIC void
vrna_ud_add_motif(vrna_fold_compound_t  *vc,
                  const char            *motif,
                  double                motif_en,
                  unsigned int          loop_type){

  if(vc){
    if(!vc->domains_up){
      /* set all default callbacks */
      vrna_ud_set_prod_rule_cb(vc, &default_prod_rule, &default_energy);
      vrna_ud_set_exp_prod_rule_cb(vc, &default_exp_prod_rule, &default_exp_energy);
      vrna_ud_set_data(vc, get_default_data(), &remove_default_data);
      vrna_ud_set_prob_cb(vc, &default_probs_add, &default_probs_get);
    }
    add_ligand_motif(vc, motif, motif_en, loop_type);
  }
}


PUBLIC int *
vrna_ud_get_motif_size_at(vrna_fold_compound_t  *vc,
                          int                   i,
                          unsigned int          loop_type){

  if(vc && vc->domains_up){
    int k, l, cnt, *ret, *ptr;

    ret = NULL;
    if((i > 0) && (i <= vc->length)){
      ptr = get_motifs(vc, i, loop_type);
      if(ptr){
        for(k = 0; ptr[k] != -1; k++) /* replace motif number with its size */
          ptr[k] = vc->domains_up->motif_size[ptr[k]];
        /* make the list unique */
        ret     = (int *)vrna_alloc(sizeof(int) * (k + 1));
        ret[0]  = -1;
        cnt     = 0;
        for(k = 0; ptr[k] != -1; k++){
          for(l = 0; l < cnt; l++){
            if(ptr[k] == ret[l])
              break; /* we've already seen this size */
          }
          if(l == cnt){ /* we've not seen this size before */
            ret[cnt]      = ptr[k];
            ret[cnt + 1]  = -1;
            cnt++;
          }
        }
        /* resize ret array */
        ret = (int *)vrna_realloc(ret, sizeof(int) * (cnt + 1));
      }
      free(ptr);
    }
    return ret;
  }
  return NULL;
}

PUBLIC int *
vrna_ud_get_motifs_at(vrna_fold_compound_t  *vc,
                      int                   i,
                      unsigned int          loop_type){

  if(vc && vc->domains_up){
    if((i > 0) && (i <= vc->length)){
      return get_motifs(vc, i, loop_type);
    }
  }
  return NULL;
}

vrna_ud_motif_t *
vrna_ud_detect_motifs(vrna_fold_compound_t  *vc,
                      const char            *structure){

  int list_size, list_pos;
  vrna_ud_motif_t *motif_list;

  motif_list = NULL;

  if(structure && vc->domains_up){
    int   l, start, end;
    char  last, *loops;

    l           = 0;
    list_pos    = 0;
    list_size   = 15;
    motif_list  = (vrna_ud_motif_t *)vrna_alloc(sizeof(vrna_ud_motif_t) * list_size);
    loops       = vrna_db_to_element_string(structure);

    while(l < vc->length){
      /* skip uppercase encodings */
      while(l < vc->length){
        if(islower(loops[l]))
          break;
        l++;
      }

      if(l < vc->length){
        start = 1 + l;
        last  = loops[l];
        while(loops[l++] == last){
          if(l == vc->length)
            break;
        }
        end = l - 1;
        annotate_ud(vc, start, end, last, &motif_list, &list_size, &list_pos);
      }
    }

    motif_list = (vrna_ud_motif_t *)vrna_realloc(motif_list, sizeof(vrna_ud_motif_t) * (list_pos + 1));
    motif_list[list_pos].start = 0;
    motif_list[list_pos].number = -1;
    free(loops);
  }

  return motif_list;
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
  data->exp_dG            = NULL;
  data->energies_ext      = NULL;
  data->energies_hp       = NULL;
  data->energies_int      = NULL;
  data->energies_mb       = NULL;
  data->exp_energies_ext  = NULL;
  data->exp_energies_hp   = NULL;
  data->exp_energies_int  = NULL;
  data->exp_energies_mb   = NULL;
  data->outside_ext       = NULL;
  data->outside_hp        = NULL;
  data->outside_int       = NULL;
  data->outside_mb        = NULL;
  data->outside_ext_count = NULL;
  data->outside_hp_count  = NULL;
  data->outside_int_count = NULL;
  data->outside_mb_count  = NULL;
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

  free(vc->domains_up->uniq_motif_size);

  free(vc->domains_up);

  vc->domains_up = NULL;
}

PRIVATE void
init_ligands_up(vrna_fold_compound_t *vc){

  vc->domains_up = (vrna_ud_t *)vrna_alloc(sizeof(vrna_ud_t));

  vc->domains_up->uniq_motif_count  = 0;
  vc->domains_up->uniq_motif_size   = NULL;
  vc->domains_up->motif_count       = 0;
  vc->domains_up->motif             = NULL;
  vc->domains_up->motif_size        = NULL;
  vc->domains_up->motif_en          = NULL;
  vc->domains_up->motif_type        = NULL;
  vc->domains_up->prod_cb           = NULL;
  vc->domains_up->exp_prod_cb       = NULL;
  vc->domains_up->energy_cb         = NULL;
  vc->domains_up->exp_energy_cb     = NULL;
  vc->domains_up->data              = NULL;
  vc->domains_up->free_data         = NULL;
  vc->domains_up->probs_add         = NULL;
  vc->domains_up->probs_get         = NULL;
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

  unsigned int  i, n, same_size;
  vrna_ud_t     *ud;

  n   = (unsigned int)strlen(motif);
  ud  = vc->domains_up;

  /* First, we update the list of unique motif lengths */
  for(same_size = i = 0; i < ud->uniq_motif_count; i++){
    if(ud->uniq_motif_size[i] == n){
      same_size = 1;
      break;
    }
  }

  if(!same_size){
    ud->uniq_motif_size = (unsigned int *)vrna_realloc(ud->uniq_motif_size, sizeof(unsigned int *) * (ud->uniq_motif_count + 1));
    ud->uniq_motif_size[ud->uniq_motif_count] = n;
    ud->uniq_motif_count++;
  }

  /* And finally, we store the motif */
  ud->motif                       = (char **)vrna_realloc(ud->motif, sizeof(char *) * (ud->motif_count + 1));
  ud->motif[ud->motif_count]      = strdup(motif);
  ud->motif_size                  = (unsigned int *)vrna_realloc(ud->motif_size, sizeof(unsigned int *) * (ud->motif_count + 1));
  ud->motif_size[ud->motif_count] = n;
  ud->motif_en                    = (double *)vrna_realloc(ud->motif_en, sizeof(double) * (ud->motif_count + 1));
  ud->motif_en[ud->motif_count]   = motif_en;
  ud->motif_type                  = (unsigned int *)vrna_realloc(ud->motif_type, sizeof(double) * (ud->motif_count + 1));
  ud->motif_type[ud->motif_count] = loop_type;
  ud->motif_count++;

}

PRIVATE void
remove_default_data(void *d){

  struct ligands_up_data_default  *data;

  data = (struct ligands_up_data_default *)d;

  free_default_data_matrices(data);
  free_default_data_exp_matrices(data);
  free_default_data(data);

  free(data->dG);
  free(data->exp_dG);
}

PRIVATE void
free_default_data(struct ligands_up_data_default *data){

  int i;
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
}

PRIVATE void
free_default_data_exp_matrices(struct ligands_up_data_default *data){

  int i;

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
    data->exp_energies_ext = NULL;
  }
  if(data->exp_energies_hp){
    /* check whether one of the other b* points to the same memory location */
    if(data->exp_energies_hp == data->exp_energies_int)
      data->exp_energies_int = NULL;
    if(data->exp_energies_hp == data->exp_energies_mb)
      data->exp_energies_mb = NULL;
    free(data->exp_energies_hp);
    data->exp_energies_hp = NULL;
  }
  if(data->exp_energies_int){
    /* check whether one of the other b* points to the same memory location */
    if(data->exp_energies_int == data->exp_energies_mb)
      data->exp_energies_mb = NULL;
    free(data->exp_energies_int);
    data->exp_energies_int = NULL;
  }
  free(data->exp_energies_mb);
  data->exp_energies_mb = NULL;

  if(data->outside_ext)
    for(i = 0; i <= data->n; i++)
      if(data->outside_ext[i])
        free(data->outside_ext[i]);
  free(data->outside_ext);
  free(data->outside_ext_count);

  if(data->outside_hp)
    for(i = 0; i <= data->n; i++)
      if(data->outside_hp[i])
        free(data->outside_hp[i]);
  free(data->outside_hp);
  free(data->outside_hp_count);

  if(data->outside_int)
    for(i = 0; i <= data->n; i++)
      if(data->outside_int[i])
        free(data->outside_int[i]);
  free(data->outside_int);
  free(data->outside_int_count);

  if(data->outside_mb)
    for(i = 0; i <= data->n; i++)
      if(data->outside_mb[i])
        free(data->outside_mb[i]);
  free(data->outside_mb);
  free(data->outside_mb_count);
}

PRIVATE int *
get_motifs(vrna_fold_compound_t *vc, int i, unsigned int loop_type){

  int               k, j, u, n, *motif_list, cnt, guess;
  char              *sequence;
  vrna_ud_t         *domains_up;

  sequence    = vc->sequence;
  n           = (int)vc->length;
  domains_up  = vc->domains_up;

  cnt         = 0;
  guess       = domains_up->motif_count;
  motif_list  = (int *)vrna_alloc(sizeof(int) * (guess + 1));

  /* collect list of motif numbers we find that start at position i */
  for(k = 0; k < domains_up->motif_count; k++){

    if(!(domains_up->motif_type[k] & loop_type))
      continue;

    j = i + domains_up->motif_size[k] - 1;
    if(j <= n){ /* only consider motif that does not exceed sequence length (does not work for circular RNAs!) */
      for(u = i; u <= j; u++){
        if(!vrna_nucleotide_IUPAC_identity(sequence[u-1], domains_up->motif[k][u-i]))
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


static void
annotate_ud(vrna_fold_compound_t *vc,
            int start,
            int end,
            char l,
            vrna_ud_motif_t **list,
            int *list_size,
            int *list_pos){

  int i,j;

  /* get motifs in segment [start,end] */
  for(i = start; i <= end; i++){
    unsigned int type = 0;
    switch(l){
      case 'e': type = VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP;
                break;
      case 'h': type = VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP;
                break;
      case 'i': type = VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP;
                break;
      case 'm': type = VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP;
                break;
    }

    int *m = vrna_ud_get_motifs_at(vc, i, type);
    if(m){
      for(j = 0; m[j] != -1; j++){
        int size = vc->domains_up->motif_size[m[j]];

        if(i + size - 1 <= end){
          if(*list_pos == *list_size){
            *list_size *= 1.2;
            *list = (vrna_ud_motif_t *)vrna_realloc(*list, sizeof(vrna_ud_motif_t) * (*list_size));
          }
          (*list)[*list_pos].start = i;
          (*list)[*list_pos].number = m[j];
          (*list_pos)++;
        }
      }
    }
    free(m);
  }
}


PRIVATE void
prepare_matrices( vrna_fold_compound_t *vc,
                  struct ligands_up_data_default *data){

  int         i,j,k,n,size;
  vrna_ud_t   *domains_up;

  n             = (int)vc->length;
  size          = ((n+1)*(n+2))/2 + 1;
  domains_up    = vc->domains_up;

  free_default_data_matrices(data);

  /* here we save memory by re-using DP matrices */
  unsigned int lt[4] = {  VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                          VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                          VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                          VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP };
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
    col     = (unsigned int *)vrna_alloc(sizeof(unsigned int) * domains_up->motif_count);
    col2    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * domains_up->motif_count);
    *(m[i]) = mx;

    for(k = 0; k < domains_up->motif_count; k++)
      col[k] = domains_up->motif_type[k] & lt[i];

    /* check if any of the remaining DP matrices can point to the same location */
    for(j=i+1;j<4;j++){
      for(k = 0; k < domains_up->motif_count; k++){
        col2[k] = domains_up->motif_type[k] & lt[j];
        if(col[k] != col2[k])
          break;
      }
      if(k == domains_up->motif_count){
        *(m[j]) = mx;
      }
    }

    free(col);
    free(col2);
  }
}

PRIVATE void
prepare_exp_matrices( vrna_fold_compound_t *vc,
                      struct ligands_up_data_default *data){

  int         i,j,k,n,size;
  vrna_ud_t   *domains_up;

  n             = (int)vc->length;
  size          = ((n+1)*(n+2))/2 + 1;
  domains_up    = vc->domains_up;

  free_default_data_exp_matrices(data);

  /* here we save memory by re-using DP matrices */
  unsigned int lt[4] = {  VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP,
                          VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                          VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
                          VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP };
  FLT_OR_DBL **m[4], *mx;
  m[0] = &data->exp_energies_ext;
  m[1] = &data->exp_energies_hp;
  m[2] = &data->exp_energies_int;
  m[3] = &data->exp_energies_mb;

  for(i=0; i<4; i++){
    unsigned int *col,*col2;
    if(*(m[i]))
      continue;

    mx      = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * size);
    col     = (unsigned int *)vrna_alloc(sizeof(unsigned int) * domains_up->motif_count);
    col2    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * domains_up->motif_count);
    *(m[i]) = mx;

    for(k = 0; k < domains_up->motif_count; k++)
      col[k] = domains_up->motif_type[k] & lt[i];

    /* check if any of the remaining DP matrices can point to the same location */
    for(j=i+1;j<4;j++){
      for(k = 0; k < domains_up->motif_count; k++){
        col2[k] = domains_up->motif_type[k] & lt[j];
        if(col[k] != col2[k])
          break;
      }
      if(k == domains_up->motif_count){
        *(m[j]) = mx;
      }
    }

    free(col);
    free(col2);
  }

  /* now prepate memory for outside partition function */
  data->outside_ext = (struct default_outside **)vrna_alloc(sizeof(struct default_outside *) * (n + 2));
  data->outside_hp  = (struct default_outside **)vrna_alloc(sizeof(struct default_outside *) * (n + 2));
  data->outside_int = (struct default_outside **)vrna_alloc(sizeof(struct default_outside *) * (n + 2));
  data->outside_mb  = (struct default_outside **)vrna_alloc(sizeof(struct default_outside *) * (n + 2));
  data->outside_ext_count = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
  data->outside_hp_count  = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
  data->outside_int_count = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
  data->outside_mb_count  = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
}

PRIVATE void
prepare_default_data( vrna_fold_compound_t *vc,
                      struct ligands_up_data_default *data){

  int         i, n;
  vrna_ud_t   *domains_up;

  n           = (int)vc->length;
  domains_up  = vc->domains_up;

  data->n = n;
  free_default_data(data);

  /*
      create motif_list for associating a nucleotide position with all
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
    data->motif_list_mb[i]  = get_motifs(vc, i, VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP);
  }

  /*  store length of motifs in 'data' */
  data->len = (int *)vrna_alloc(sizeof(int) * domains_up->motif_count);
  for(i = 0; i < domains_up->motif_count; i++)
    data->len[i] = domains_up->motif_size[i];

}

PRIVATE void
default_prod_rule(vrna_fold_compound_t *vc,
                  void *d){

  int                             i,j,k,l,u,n,size,e_ext, e_hp, e_int, e_mb,en,en2,*idx;
  unsigned int                    loop_type;
  vrna_ud_t                       *domains_up;
  struct ligands_up_data_default  *data;

  int           *energies_ext;
  int           *energies_hp;
  int           *energies_int;
  int           *energies_mb;

  n             = (int)vc->length;
  size          = ((n+1)*(n+2))/2 + 1;
  idx           = vc->jindx;
  domains_up    = vc->domains_up;
  data          = (struct ligands_up_data_default *)d;

  prepare_default_data(vc, data);
  prepare_matrices(vc, data);

  energies_ext  = data->energies_ext;
  energies_hp   = data->energies_hp;
  energies_int  = data->energies_int;
  energies_mb   = data->energies_mb;

  data->default_cb[VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP] = default_exp_energy_ext_motif;
  data->default_cb[VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP]  = default_exp_energy_hp_motif;
  data->default_cb[VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP] = default_exp_energy_int_motif;
  data->default_cb[VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP]  = default_exp_energy_mb_motif;

  /*  precompute energy contributions of the motifs */
  data->dG = (int *)vrna_alloc(sizeof(int) * domains_up->motif_count);
  for(i = 0; i < domains_up->motif_count; i++)
    data->dG[i] = (int)roundf(domains_up->motif_en[i] * 100.);

  /* now we can start to fill the DP matrices */
  for(i=n; i>0; i--){
    int *list_ext = data->motif_list_ext[i];
    int *list_hp  = data->motif_list_hp[i];
    int *list_int = data->motif_list_int[i];
    int *list_mb  = data->motif_list_mb[i];
    for(j=i;j<=n;j++){
      if(i < j){
        e_ext = energies_ext[idx[j]+i+1];
        e_hp  = energies_hp[idx[j]+i+1];
        e_int = energies_int[idx[j]+i+1];
        e_mb  = energies_mb[idx[j]+i+1];
      } else {
        e_ext = INF;
        e_hp  = INF;
        e_int = INF;
        e_mb  = INF;
      }
      if(list_ext){
        for(k = 0; -1 != (l = list_ext[k]); k++){
          u   = i + data->len[l] - 1;
          en  = data->dG[l];
          if(u <= j){
            e_ext = MIN2(e_ext, en);
            if(u < j){
              en2 = en + energies_ext[idx[j]+u+1];
              e_ext = MIN2(e_ext, en2);
            }
          }
        }
      }
      if(list_hp){
        for(k = 0; -1 != (l = list_hp[k]); k++){
          u   = i + data->len[l] - 1;
          en  = data->dG[l];
          if(u <= j){
            e_hp = MIN2(e_hp, en);
            if(u < j){
              en2 = en + energies_hp[idx[j]+u+1];
              e_hp = MIN2(e_hp, en2);
            }
          }
        }
      }
      if(list_int){
        for(k = 0; -1 != (l = list_int[k]); k++){
          u   = i + data->len[l] - 1;
          en  = data->dG[l];
          if(u <= j){
            e_int = MIN2(e_int, en);
            if(u < j){
              en2 = en + energies_int[idx[j]+u+1];
              e_int = MIN2(e_int, en2);
            }
          }
        }
      }
      if(list_mb){
        for(k = 0; -1 != (l = list_mb[k]); k++){
          u   = i + data->len[l] - 1;
          en  = data->dG[l];
          if(u <= j){
            e_mb = MIN2(e_mb, en);
            if( u < j){
              en2 = en + energies_mb[idx[j]+u+1];
              e_mb = MIN2(e_mb, en2);
            }
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

  int                             i,j,k,l,u,n,size,*idx;
  unsigned int                    loop_type;
  FLT_OR_DBL                      q_ext, q_hp, q_int, q_mb, q, qq;
  vrna_ud_t                       *domains_up;
  struct ligands_up_data_default  *data;

  FLT_OR_DBL    *exp_energies_ext;
  FLT_OR_DBL    *exp_energies_hp;
  FLT_OR_DBL    *exp_energies_int;
  FLT_OR_DBL    *exp_energies_mb;
  double        kT;

  n             = (int)vc->length;
  size          = ((n+1)*(n+2))/2 + 1;
  idx           = vc->iindx;
  domains_up    = vc->domains_up;
  data          = (struct ligands_up_data_default *)d;
  kT            = vc->exp_params->kT;

  prepare_default_data(vc, data);
  prepare_exp_matrices(vc, data);

  exp_energies_ext  = data->exp_energies_ext;
  exp_energies_hp   = data->exp_energies_hp;
  exp_energies_int  = data->exp_energies_int;
  exp_energies_mb   = data->exp_energies_mb;

  data->default_cb[VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP] = default_exp_energy_ext_motif;
  data->default_cb[VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP]  = default_exp_energy_hp_motif;
  data->default_cb[VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP] = default_exp_energy_int_motif;
  data->default_cb[VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP]  = default_exp_energy_mb_motif;

  data->exp_e_mx[VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP] = data->exp_energies_ext;
  data->exp_e_mx[VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP]  = data->exp_energies_hp;
  data->exp_e_mx[VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP] = data->exp_energies_int;
  data->exp_e_mx[VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP]  = data->exp_energies_mb;

  /*  precompute energy contributions of the motifs */
  data->exp_dG = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * domains_up->motif_count);
  for(i = 0; i < domains_up->motif_count; i++){
    double GT = domains_up->motif_en[i] * 1000.; /* in cal/mol */
    data->exp_dG[i] = (FLT_OR_DBL)exp( -GT / kT);
  }

  /* now we can start to fill the DP matrices */
  for(i=n; i>0; i--){
    int *list_ext = data->motif_list_ext[i];
    int *list_hp  = data->motif_list_hp[i];
    int *list_int = data->motif_list_int[i];
    int *list_mb  = data->motif_list_mb[i];
    for(j = i; j <= n; j++){
      if(i < j){
        q_ext = exp_energies_ext[idx[i + 1] - j];
        q_hp  = exp_energies_hp[idx[i + 1] - j];
        q_int = exp_energies_int[idx[i + 1] - j];
        q_mb  = exp_energies_mb[idx[i + 1] - j];
      } else {
        q_ext = 0;
        q_hp  = 0;
        q_int = 0;
        q_mb  = 0;
      }
      if(list_ext){
        for(k = 0; -1 != (l = list_ext[k]); k++){
          u = i + data->len[l] - 1;
          q = data->exp_dG[l];
          if(u <= j){
            q_ext += q;
            if(u < j)
              q_ext += q * exp_energies_ext[idx[u + 1] - j];
          }
        }
      }
      if(list_hp){
        for(k = 0; -1 != (l = list_hp[k]); k++){
          u   = i + data->len[l] - 1;
          q  = data->exp_dG[l];
          if(u <= j){
            q_hp += q;
            if(u < j)
              q_hp += q * exp_energies_hp[idx[u + 1] - j];
          }
        }
      }
      if(list_int){
        for(k = 0; -1 != (l = list_int[k]); k++){
          u   = i + data->len[l] - 1;
          q  = data->exp_dG[l];
          if(u <= j){
            q_int += q;
            if(u < j)
              q_int += q * exp_energies_int[idx[u + 1] - j];
          }
        }
      }
      if(list_mb){
        for(k = 0; -1 != (l = list_mb[k]); k++){
          u   = i + data->len[l] - 1;
          q  = data->exp_dG[l];
          if(u <= j){
            q_mb += q;
            if(u < j)
              q_mb += q * exp_energies_mb[idx[u + 1] - j];
          }
        }
      }

      exp_energies_ext[idx[i] - j]  = q_ext;
      exp_energies_hp[idx[i] - j]   = q_hp;
      exp_energies_int[idx[i] - j]  = q_int;
      exp_energies_mb[idx[i] - j]   = q_mb;
    }
  }
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
    return INF;

  if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_MOTIF){
    if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP)
      en = default_energy_ext_motif(i, j, data);
    else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP)
      en = default_energy_hp_motif(i, j, data);
    else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP)
      en = default_energy_int_motif(i, j, data);
    else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP)
      en = default_energy_mb_motif(i, j, data);
  } else {
    if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP){
      if(data->energies_ext)
        en = data->energies_ext[ij];
    } else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP){
      if(data->energies_hp)
        en = data->energies_hp[ij];
    } else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP){
      if(data->energies_int)
        en = data->energies_int[ij];
    } else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP){
      if(data->energies_mb)
        en = data->energies_mb[ij];
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

  if(j < i)
    return 0.;

  if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_MOTIF){
    q = data->default_cb[loop_type & ~(VRNA_UNSTRUCTURED_DOMAIN_MOTIF)](i, j, data);
  } else {
    idx = vc->iindx;
    ij  = idx[i] - j;
    q   = data->exp_e_mx[loop_type][ij];
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
        e = MIN2(e, data->dG[m]);
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
        e = MIN2(e, data->dG[m]);
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
        e = MIN2(e, data->dG[m]);
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
        e = MIN2(2, data->dG[m]);
      k++;
    }
  }
  return e;
}

PRIVATE FLT_OR_DBL
default_exp_energy_ext_motif( int i,
                              int j,
                              struct ligands_up_data_default *data){

  int         k, m;
  FLT_OR_DBL  q = 0;

  if(data->motif_list_ext[i]){
    k = 0;
    while(-1 != (m = data->motif_list_ext[i][k])){
      if((i + data->len[m] - 1) == j)
        q += data->exp_dG[m];
      k++;
    }
  }

  return q;
}

PRIVATE FLT_OR_DBL
default_exp_energy_hp_motif(int i,
                            int j,
                            struct ligands_up_data_default *data){

  int         k, m;
  FLT_OR_DBL  q = 0;

  if(data->motif_list_hp[i]){
    k = 0;
    while(-1 != (m = data->motif_list_hp[i][k])){
      if((i + data->len[m] - 1) == j)
        q += data->exp_dG[m];
      k++;
    }
  }

  return q;
}

PRIVATE FLT_OR_DBL
default_exp_energy_int_motif( int i,
                              int j,
                              struct ligands_up_data_default *data){

  int         k, m;
  FLT_OR_DBL  q = 0;

  if(data->motif_list_int[i]){
    k = 0;
    while(-1 != (m = data->motif_list_int[i][k])){
      if((i + data->len[m] - 1) == j)
        q += data->exp_dG[m];
      k++;
    }
  }

  return q;
}

PRIVATE FLT_OR_DBL
default_exp_energy_mb_motif(int i,
                            int j,
                            struct ligands_up_data_default *data){

  int         k, m;
  FLT_OR_DBL  q = 0;

  if(data->motif_list_mb[i]){
    k = 0;
    while(-1 != (m = data->motif_list_mb[i][k])){
      if((i + data->len[m] - 1) == j)
        q += data->exp_dG[m];
      k++;
    }
  }

  return q;
}

PRIVATE void
default_probs_add(vrna_fold_compound_t *vc,
                  int i,
                  int j,
                  unsigned int loop_type,
                  FLT_OR_DBL exp_energy,
                  void *data){

  int                             **motif_list, k, l, m;
  unsigned int                    *size, *cnt, o;
  struct ligands_up_data_default  *d;
  struct default_outside          **storage, **st;

  d = (struct ligands_up_data_default  *)data;

  if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_MOTIF){
    if(j < i)
      return;

    if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP){
      motif_list  = d->motif_list_ext;
      storage     = &(d->outside_ext[i]);
      size        = &(d->outside_ext_count[i]);
    } else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP) {
      motif_list = d->motif_list_hp;
      storage     = &(d->outside_hp[i]);
      size        = &(d->outside_hp_count[i]);
    } else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP){
      motif_list = d->motif_list_int;
      storage     = &(d->outside_int[i]);
      size        = &(d->outside_int_count[i]);
    } else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP){
      motif_list = d->motif_list_mb;
      storage     = &(d->outside_mb[i]);
      size        = &(d->outside_mb_count[i]);
    } else{
      vrna_message_warning("Unknown unstructured domain loop type");
      return;
    }

    k = 0;
    while(-1 != (m = motif_list[i][k])){
      if((i + d->len[m] - 1) == j){

        /* check for addition first */
        for(o = 0; o < *size; o++)
          if((*storage)[o].motif_num == m){ /* found previously added motif constribution */
            (*storage)[o].exp_energy += exp_energy;
            break;
          }

        /* if we haven't added yet, create new list entry */
        if(o == *size){
          *storage = (struct default_outside *)vrna_realloc(*storage, sizeof(struct default_outside) * (*size + 1));
          (*storage)[*size].motif_num   = m;
          (*storage)[*size].exp_energy  = exp_energy;
          (*size)++;
        }
      }
      k++;
    }

  } else {
    if(j < i)
      return;

    FLT_OR_DBL pf, exp_e;
    pf = default_exp_energy(vc, i, j, loop_type, data);

    if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP){
      motif_list  = d->motif_list_ext;
      storage     = d->outside_ext;
      size        = d->outside_ext_count;
    } else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP) {
      motif_list = d->motif_list_hp;
      storage     = d->outside_hp;
      size        = d->outside_hp_count;
    } else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP){
      motif_list = d->motif_list_int;
      storage     = d->outside_int;
      size        = d->outside_int_count;
    } else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP){
      motif_list = d->motif_list_mb;
      storage     = d->outside_mb;
      size        = d->outside_mb_count;
    } else{
      vrna_message_warning("Unknown unstructured domain loop type");
      return;
    }

    /* check for each motif starting at any k with i <= k <= j */
    for(k = i; k <= j; k++){
      if(motif_list[k]){
        st  = &(storage[k]);
        cnt = &(size[k]);
        for(l = 0; motif_list[k][l] != -1; l++){
          m = motif_list[k][l];

          if(j < k + d->len[m] - 1) /* motifs must be sorted be length */
            continue;

          exp_e         = d->exp_dG[m];
          FLT_OR_DBL p  = exp_e / pf;

          /* add/insert contribution */

          /* check for addition first */
          for(o = 0; o < *cnt; o++)
            if((*st)[o].motif_num == m){ /* found previously added motif constribution */
              (*st)[o].exp_energy += p * exp_energy;
              break;
            }

          /* if we haven't added yet, create new list entry */
          if(o == *cnt){
            *st = (struct default_outside *)vrna_realloc(*st, sizeof(struct default_outside) * (*cnt + 1));
            (*st)[*cnt].motif_num   = m;
            (*st)[*cnt].exp_energy  = p * exp_energy;
            (*cnt)++;
          }
        }
      }
    }
  }

}

PRIVATE FLT_OR_DBL
default_probs_get(vrna_fold_compound_t *vc,
                  int i,
                  int j,
                  unsigned int loop_type,
                  int motif,
                  void *data){

  FLT_OR_DBL                      outside = 0.;
  unsigned int                    *size, k;
  struct ligands_up_data_default  *d;
  struct default_outside          **storage;

  d = (struct ligands_up_data_default  *)data;

  if(j < i)
    return 0.;

  if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP){
    storage     = &(d->outside_ext[i]);
    size        = &(d->outside_ext_count[i]);
  } else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP) {
    storage     = &(d->outside_hp[i]);
    size        = &(d->outside_hp_count[i]);
  } else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP){
    storage     = &(d->outside_int[i]);
    size        = &(d->outside_int_count[i]);
  } else if(loop_type & VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP){
    storage     = &(d->outside_mb[i]);
    size        = &(d->outside_mb_count[i]);
  } else{
    vrna_message_warning("Unknown unstructured domain loop type");
    return 0.;
  }

  for(k = 0; k < *size; k++){
    /* check for motif number match */
    if((*storage)[k].motif_num == motif)
      /* check for length match */
      if(i + d->len[motif] - 1 == j){
        outside += (*storage)[k].exp_energy;
      }
  }

  return outside;
}
