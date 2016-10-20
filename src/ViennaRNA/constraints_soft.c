/* constraints handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/energy_const.h" /* defines MINPSCORE */
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/aln_util.h"
#include "ViennaRNA/file_formats.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/constraints_soft.h"


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

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE void
sc_add_up(vrna_fold_compound_t *vc,
          const FLT_OR_DBL *constraints);

PRIVATE void
sc_add_bp(vrna_fold_compound_t *vc,
          const FLT_OR_DBL **constraints);

PRIVATE void
sc_really_add_up( vrna_fold_compound_t *vc,
                  int i,
                  FLT_OR_DBL energy);

PRIVATE void
sc_really_add_bp( vrna_fold_compound_t *vc,
                  int i,
                  int j,
                  FLT_OR_DBL energy);

PRIVATE void
prepare_Boltzmann_weights_up( vrna_fold_compound_t *vc);

PRIVATE void
prepare_Boltzmann_weights_bp(vrna_fold_compound_t *vc);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC void
vrna_sc_init(vrna_fold_compound_t *vc){

  unsigned int s;
  vrna_sc_t    *sc;

  if(vc){
    vrna_sc_remove(vc);

    switch(vc->type){
      case VRNA_FC_TYPE_SINGLE:     sc                    = (vrna_sc_t *)vrna_alloc(sizeof(vrna_sc_t));
                                    sc->energy_up         = NULL;
                                    sc->energy_bp         = NULL;
                                    sc->energy_stack      = NULL;
                                    sc->exp_energy_stack  = NULL;
                                    sc->exp_energy_up     = NULL;
                                    sc->exp_energy_bp     = NULL;
                                    sc->f                 = NULL;
                                    sc->exp_f             = NULL;
                                    sc->data              = NULL;
                                    sc->free_data         = NULL;

                                    vc->sc  = sc;
                                    break;

      case VRNA_FC_TYPE_COMPARATIVE:  vc->scs = (vrna_sc_t **)vrna_alloc(sizeof(vrna_sc_t*) * (vc->n_seq + 1));
                                    for(s = 0; s < vc->n_seq; s++){
                                      sc                    = (vrna_sc_t *)vrna_alloc(sizeof(vrna_sc_t));
                                      sc->energy_up         = NULL;
                                      sc->energy_bp         = NULL;
                                      sc->energy_stack      = NULL;
                                      sc->exp_energy_stack  = NULL;
                                      sc->exp_energy_up     = NULL;
                                      sc->exp_energy_bp     = NULL;
                                      sc->f                 = NULL;
                                      sc->exp_f             = NULL;
                                      sc->data              = NULL;
                                      sc->free_data         = NULL;

                                      vc->scs[s]  = sc;
                                    }
                                    break;
      default:                      /* do nothing */
                                    break;
    }
  }
}

PUBLIC void
vrna_sc_remove(vrna_fold_compound_t *vc){

  int s;

  if(vc){
    switch(vc->type){
      case  VRNA_FC_TYPE_SINGLE:    vrna_sc_free(vc->sc);
                                    vc->sc = NULL;
                                    break;
      case  VRNA_FC_TYPE_COMPARATIVE: if(vc->scs){
                                      for(s = 0; s < vc->n_seq; s++)
                                        vrna_sc_free(vc->scs[s]);
                                      free(vc->scs);
                                    }
                                    vc->scs = NULL;
                                    break;
      default:                      /* do nothing */
                                    break;
    }
  }
}

PUBLIC void
vrna_sc_free(vrna_sc_t *sc){

  int i;
  if(sc){
    if(sc->energy_up){
      for(i = 0; sc->energy_up[i]; free(sc->energy_up[i++]));
      free(sc->energy_up);
    }
    if(sc->exp_energy_up){
      for(i = 0; sc->exp_energy_up[i]; free(sc->exp_energy_up[i++]));
      free(sc->exp_energy_up);
    }
    
    free(sc->energy_bp);
    free(sc->exp_energy_bp);
    free(sc->energy_stack);
    free(sc->exp_energy_stack);

    if(sc->free_data)
      sc->free_data(sc->data);

    free(sc);
  }
}

PUBLIC void
vrna_sc_set_bp( vrna_fold_compound_t *vc,
                const FLT_OR_DBL **constraints,
                unsigned int options){

  if(vc && (vc->type == VRNA_FC_TYPE_SINGLE)){
    if(constraints){
      /* always add (pure) soft constraints */
      sc_add_bp(vc, constraints);
    }

    if(options & VRNA_OPTION_PF) /* prepare Boltzmann factors for the BP soft constraints */
      prepare_Boltzmann_weights_bp(vc);
  }
}

PUBLIC void
vrna_sc_add_bp(vrna_fold_compound_t *vc,
                      int i,
                      int j,
                      FLT_OR_DBL energy,
                      unsigned int options){

  if(vc && (vc->type == VRNA_FC_TYPE_SINGLE)){
    sc_really_add_bp(vc, i, j, energy);

    if(options & VRNA_OPTION_PF) /* prepare Boltzmann factors for the BP soft constraints */
      prepare_Boltzmann_weights_bp(vc);
  }
}

PUBLIC void
vrna_sc_set_up( vrna_fold_compound_t *vc,
                const FLT_OR_DBL *constraints,
                unsigned int options){

  if(vc && (vc->type == VRNA_FC_TYPE_SINGLE)){
    if(constraints){
      /* always add (pure) soft constraints */
      sc_add_up(vc, constraints);
    }

    if(options & VRNA_OPTION_PF)
      prepare_Boltzmann_weights_up(vc);
  }
}

PUBLIC void
vrna_sc_add_up(vrna_fold_compound_t *vc,
                      int i,
                      FLT_OR_DBL energy,
                      unsigned int options){

  if(vc && (vc->type == VRNA_FC_TYPE_SINGLE)){
    if((i < 1) || (i > vc->length)){
      vrna_message_warning( "vrna_sc_add_up(): Nucleotide position %d out of range!"
                            " (Sequence length: %d)",
                            i, vc->length);
    } else {
      sc_really_add_up(vc, i, energy);

      if(options & VRNA_OPTION_PF)
        prepare_Boltzmann_weights_up(vc);
    }
  }
}

PUBLIC void
vrna_sc_add_data( vrna_fold_compound_t *vc,
                  void *data,
                  vrna_callback_free_auxdata *free_data){

  if(vc){
    if(vc->type == VRNA_FC_TYPE_SINGLE){
      if(!vc->sc)
        vrna_sc_init(vc);

      vc->sc->data        = data;
      vc->sc->free_data   = free_data;
    }
  }
}

PUBLIC void
vrna_sc_add_f(vrna_fold_compound_t *vc,
              vrna_callback_sc_energy *f){

  if(vc && f){
    if(vc->type == VRNA_FC_TYPE_SINGLE){
      if(!vc->sc)
        vrna_sc_init(vc);

      vc->sc->f       = f;
    }
  }
}

PUBLIC void
vrna_sc_add_bt( vrna_fold_compound_t *vc,
                vrna_callback_sc_backtrack *f){

  if(vc && f){
    if(vc->type == VRNA_FC_TYPE_SINGLE){
      if(!vc->sc)
        vrna_sc_init(vc);

      vc->sc->bt      = f;
    }
  }
}

PUBLIC void
vrna_sc_add_exp_f(vrna_fold_compound_t *vc,
                  vrna_callback_sc_exp_energy *exp_f){

  if(vc && exp_f){
    if(vc->type == VRNA_FC_TYPE_SINGLE){
      if(!vc->sc)
        vrna_sc_init(vc);

      vc->sc->exp_f   = exp_f;
    }
  }
}

PRIVATE void
sc_add_bp(vrna_fold_compound_t *vc,
          const FLT_OR_DBL **constraints){

  unsigned int  i, j, n;
  vrna_sc_t     *sc;
  int           *idx;

  n   = vc->length;

  if(!vc->sc)
    vrna_sc_init(vc);

  sc              = vc->sc;

  if(sc->energy_bp)
    free(sc->energy_bp);

  sc->energy_bp = (int *)vrna_alloc(sizeof(int) * (((n + 1) * (n + 2)) / 2));

  idx = vc->jindx;
  for(i = 1; i < n; i++)
    for(j=i+1; j<=n; j++)
      sc->energy_bp[idx[j]+i] = (int)roundf(constraints[i][j] * 100.);

}

PRIVATE void
sc_really_add_bp( vrna_fold_compound_t *vc,
                  int i,
                  int j,
                  FLT_OR_DBL energy){

  unsigned int  n;
  vrna_sc_t     *sc;
  int           *idx;

  n   = vc->length;

  if(!vc->sc)
    vrna_sc_init(vc);

  sc              = vc->sc;

  if(!sc->energy_bp)
    sc->energy_bp = (int *)vrna_alloc(sizeof(int) * (((n + 1) * (n + 2)) / 2));

  idx = vc->jindx;
  sc->energy_bp[idx[j]+i] += (int)roundf(energy * 100.);
}


PRIVATE void
sc_add_up(vrna_fold_compound_t *vc,
          const FLT_OR_DBL *constraints){

  unsigned int  i, j, n;
  vrna_sc_t     *sc;

  if(vc && constraints){
    n   = vc->length;

    if(!vc->sc)
      vrna_sc_init(vc);

    sc  = vc->sc;
    /*  allocate memory such that we can access the soft constraint
        energies of a subsequence of length j starting at position i
        via sc->energy_up[i][j]
    */
    if(sc->energy_up){
      for(i = 0; i <= n; i++)
        if(sc->energy_up[i])
          free(sc->energy_up[i]);
      free(sc->energy_up);
    }

    sc->energy_up = (int **)vrna_alloc(sizeof(int *) * (n + 2));
    for(i = 0; i <= n; i++)
      sc->energy_up[i] = (int *)vrna_alloc(sizeof(int) * (n - i + 2));

    sc->energy_up[n+1] = NULL;

    for(i = 1; i <= n; i++){
      for(j = 1; j <= (n - i + 1); j++){
        sc->energy_up[i][j] =   sc->energy_up[i][j-1]
                                  + (int)roundf(constraints[i+j-1] * 100.); /* convert to 10kal/mol */
      }
    }
  }
}

PRIVATE void
sc_really_add_up( vrna_fold_compound_t *vc,
                  int i,
                  FLT_OR_DBL energy){

  unsigned int  j, u, n;
  vrna_sc_t     *sc;

  n   = vc->length;

  if(!vc->sc)
    vrna_sc_init(vc);

  sc  = vc->sc;
  /*
      Prepare memory:
      allocate memory such that we can access the soft constraint
      energies of a subsequence of length j starting at position i
      via sc->energy_up[i][j]
  */
  if(!(sc->energy_up)){
    sc->energy_up = (int **)vrna_alloc(sizeof(int *) * (n + 2));
  }

  for(j = 0; j <= n; j++){
    if(!(sc->energy_up[j]))
      sc->energy_up[j] = (int *)vrna_alloc(sizeof(int) * (n - j + 2));
  }

  sc->energy_up[n+1] = NULL;

  /* update soft constraints for subsequences starting at some position */
  for(j = 1; j <= i; j++){
    int u_max   = n - j + 1;
    int u_start = i - j + 1;

    sc->energy_up[j][u_start] += (int)roundf(energy * 100.); /* convert to 10kal/mol */

    for(u = u_start + 1; u <= u_max; u++){
      sc->energy_up[j][u] = sc->energy_up[j][u-1]
                            + sc->energy_up[j + u - 1][1];
    }
  }
}


/* compute Boltzmann factors for BP soft constraints from stored free energies */
PRIVATE void
prepare_Boltzmann_weights_bp(vrna_fold_compound_t *vc){

  unsigned int  i, j, n;
  vrna_sc_t     *sc;
  int           *idx, *jidx;

  if(vc->sc && vc->sc->energy_bp){
    n   = vc->length;
    sc  = vc->sc;

    vrna_exp_param_t  *exp_params = vc->exp_params;
    double            GT          = 0.;
    double            temperature = exp_params->temperature;
    double            kT          = exp_params->kT;
    double            TT          = (temperature+K0)/(Tmeasure);

    if(sc->exp_energy_bp)
      free(sc->exp_energy_bp);
    sc->exp_energy_bp     = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (((n + 1) * (n + 2)) / 2));

    idx   = vc->iindx;
    jidx  = vc->jindx;

    for(i = 1; i < n; i++)
      for(j = i + 1; j <= n; j++){
        GT = sc->energy_bp[jidx[j] + i] * 10.;
        sc->exp_energy_bp[idx[i]-j] = (FLT_OR_DBL)exp( -GT / kT);
      }
  }
}

PRIVATE void
prepare_Boltzmann_weights_up(vrna_fold_compound_t *vc){

  unsigned int  i, j, n;
  vrna_sc_t     *sc;

  n   = vc->length;

  sc  = vc->sc;

  vrna_exp_param_t   *exp_params = vc->exp_params;
  double             GT          = 0.;
  double             temperature = exp_params->temperature;
  double             kT          = exp_params->kT;
  double             TT          = (temperature+K0)/(Tmeasure);

  /* #################################### */
  /* # single nucleotide contributions  # */
  /* #################################### */

  if(sc && sc->energy_up){
    /*  allocate memory such that we can access the soft constraint
        energies of a subsequence of length j starting at position i
        via sc->exp_energy_up[i][j]
    */

    /* free previous unpaired nucleotide constraints */
    if(sc->exp_energy_up){
      for(i = 0; i <= n; i++)
        if(sc->exp_energy_up[i])
          free(sc->exp_energy_up[i]);
      free(sc->exp_energy_up);
    }

    /* allocate memory and create Boltzmann factors from stored contributions */
    sc->exp_energy_up = (FLT_OR_DBL **)vrna_alloc(sizeof(FLT_OR_DBL *) * (n + 2));
    sc->exp_energy_up[0] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n + 2));
    for(j = 0; j <= (n + 1); j++)
        sc->exp_energy_up[0][j] = 1.;

    for(i = 1; i <= n; i++){
      sc->exp_energy_up[i] = (FLT_OR_DBL *)vrna_alloc(sizeof(FLT_OR_DBL) * (n - i + 2));

      sc->exp_energy_up[i][0] = 1.;

      for(j = 1; j <= (n - i + 1); j++){
        GT  = (double)sc->energy_up[i][j] * 10.; /* convert to cal/mol */
        sc->exp_energy_up[i][j] = (FLT_OR_DBL)exp( -GT / kT);
      }
    }

    sc->exp_energy_up[n+1] = NULL;
  }
}
