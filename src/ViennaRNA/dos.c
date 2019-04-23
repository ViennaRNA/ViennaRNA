#include <stdlib.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/structures.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/datastructures/basic.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/constraints/hard.h"
#include "ViennaRNA/constraints/soft.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/structured_domains.h"
#include "ViennaRNA/unstructured_domains.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/mfe.h"

struct aux_arrays {
  int *cc;    /* auxilary arrays for canonical structures     */
  int *cc1;   /* auxilary arrays for canonical structures     */
  int *Fmi;   /* holds row i of fML (avoids jumps in memory)  */
  int *DMLi;  /* DMLi[j] holds  MIN(fML[i,k]+fML[k+1,j])      */
  int *DMLi1; /*                MIN(fML[i+1,k]+fML[k+1,j])    */
  int *DMLi2; /*                MIN(fML[i+2,k]+fML[k+1,j])    */
};


PRIVATE INLINE struct aux_arrays *
get_aux_arrays(unsigned int length)
{
  unsigned int      j;
  struct aux_arrays *aux = (struct aux_arrays *)vrna_alloc(sizeof(struct aux_arrays));

  aux->cc     = (int *)vrna_alloc(sizeof(int) * (length + 2));  /* auxilary arrays for canonical structures     */
  aux->cc1    = (int *)vrna_alloc(sizeof(int) * (length + 2));  /* auxilary arrays for canonical structures     */
  aux->Fmi    = (int *)vrna_alloc(sizeof(int) * (length + 1));  /* holds row i of fML (avoids jumps in memory)  */
  aux->DMLi   = (int *)vrna_alloc(sizeof(int) * (length + 1));  /* DMLi[j] holds  MIN(fML[i,k]+fML[k+1,j])      */
  aux->DMLi1  = (int *)vrna_alloc(sizeof(int) * (length + 1));  /*                MIN(fML[i+1,k]+fML[k+1,j])    */
  aux->DMLi2  = (int *)vrna_alloc(sizeof(int) * (length + 1));  /*                MIN(fML[i+2,k]+fML[k+1,j])    */

  /* prefill helper arrays */
  for (j = 0; j <= length; j++)
    aux->Fmi[j] = aux->DMLi[j] = aux->DMLi1[j] = aux->DMLi2[j] = INF;

  return aux;
}

PRIVATE INLINE void
rotate_aux_arrays(struct aux_arrays *aux,
                  unsigned int      length)
{
  unsigned int  j;
  int           *FF;

  FF          = aux->DMLi2;
  aux->DMLi2  = aux->DMLi1;
  aux->DMLi1  = aux->DMLi;
  aux->DMLi   = FF;
  FF          = aux->cc1;
  aux->cc1    = aux->cc;
  aux->cc     = FF;
  for (j = 1; j <= length; j++)
    aux->cc[j] = aux->Fmi[j] = aux->DMLi[j] = INF;
}


PRIVATE INLINE void
free_aux_arrays(struct aux_arrays *aux)
{
  free(aux->cc);
  free(aux->cc1);
  free(aux->Fmi);
  free(aux->DMLi);
  free(aux->DMLi1);
  free(aux->DMLi2);
  free(aux);
}


struct dp_counts_per_energy {
  int **n_ij_e;
  int **n_ij_A_e;
  int **n_ij_B_e;
  int **n_ij_M_e;
  int **n_ij_M1_e;
};

int get_energy_bin_index(int min_energy,
                         int step_energy,
                         int energy_length,
                         int target_energy){
  int bin_index = -1;
  int i = 0;
  int current_energy  = min_energy;

  int max_energy = min_energy + step_energy * energy_length;
  if(max_energy < target_energy){
    bin_index = energy_length-1; // put overflow in last bin.
  }
  else{
    for(i = 0; i < energy_length; i++){
      if(target_energy == current_energy){
        bin_index = i;
        break;
      }
      current_energy += step_energy;
    }
  }
  return bin_index;
}

PRIVATE INLINE int
decompose_pair(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               struct aux_arrays    *aux,
               int min_energy,
               int step_energy,
               int energy_length,
               struct dp_counts_per_energy *dp_count_matrix_pt)
{
  unsigned char hc_decompose;
  unsigned int  n;
  int           e, new_c, energy, stackEnergy, ij, dangle_model, noLP,
                *DMLi1, *DMLi2, *cc, *cc1;

  printf("ij %d %d \n", i,j);

  n             = fc->length;
  ij            = fc->jindx[j] + i;
  dangle_model  = fc->params->model_details.dangles;
  noLP          = fc->params->model_details.noLP;
  hc_decompose  = fc->hc->mx[n * i + j];
  DMLi1         = aux->DMLi1;
  DMLi2         = aux->DMLi2;
  cc            = aux->cc;
  cc1           = aux->cc1;
  e             = INF;
  int *rtype         = &(fc->params->model_details.rtype[0]);
  short *S1            = fc->sequence_encoding;

  int type  = fc->ptype[fc->jindx[j] + i];

  int no_close = (((type == 3) || (type == 4)) && fc->params->model_details.noGUclosure);

  /* do we evaluate this pair? */
  if (hc_decompose) {
    new_c = INF;

    /* check for hairpin loop */
    int energy_hp  = E_Hairpin(j - i - 1, type, S1[i + 1], S1[j - 1], fc->sequence + i - 1, fc->params); // vrna_E_hp_loop(fc, i, j);
    fc->matrices->c[ij] = MIN2(fc->matrices->c[ij], energy_hp);
    new_c   = MIN2(new_c, energy_hp);

    int bin_index = get_energy_bin_index(min_energy, step_energy,energy_length, energy_hp);
    printf("hp energy: %d %d\n", energy_hp, bin_index);
    dp_count_matrix_pt->n_ij_e[ij][bin_index] = 1;

    /* check for interior loops */
    //energy  = vrna_E_int_loop(fc, i, j);
    //new_c   = MIN2(new_c, energy);
    int min_int = new_c;
    /* INTERIOR LOOP STRUCTURES */
    int maxp = MIN2(j - 2 - TURN, i + MAXLOOP + 1);
    int p, q, type_2, pq;
    for (p = i + 1; p <= maxp; p++) {
      unsigned int  minq    = p + TURN + 1;

      for (q = minq; q < j; q++) {
        pq = fc->jindx[q] + p;
        /* set distance to reference structure... */
        type_2 = fc->ptype[fc->jindx[q] + p];

        if (type_2 == 0)
          continue;

        type_2 = rtype[type_2];

        if (no_closingGU)
          if (no_close || (type_2 == 3) || (type_2 == 4))
            if ((p > i + 1) || (q < j - 1))
              continue;

       /* continue unless stack */
        energy = E_IntLoop(p - i - 1, j - q - 1, type, type_2, S1[i + 1], S1[j - 1], S1[p - 1], S1[q + 1], fc->params);
        //new_c   = MIN2(new_c, energy);
        min_int = MIN2(min_int, fc->matrices->c[pq] + energy);
        if(fc->matrices->c[pq] != INF)
          fc->matrices->c[ij] = MIN2(fc->matrices->c[ij], fc->matrices->c[pq] + energy);

        int bin_index_inner_loop = get_energy_bin_index(min_energy, step_energy,energy_length,  fc->matrices->c[pq]);
        bin_index = get_energy_bin_index(min_energy, step_energy,energy_length,  energy_hp + energy);
        if(dp_count_matrix_pt->n_ij_e[pq] != NULL){
          printf("i,j,p,q, bin_index, bin_index_inner_loop, energy %d %d %d %d %d %d %d\n", i,j,p,q, bin_index, bin_index_inner_loop,  energy_hp + energy);
          printf("e_int, e_pq %d %d \n", energy, fc->matrices->c[pq]);
          if(dp_count_matrix_pt->n_ij_e[pq] != NULL){
            //printf("address from: %#08x \n", &dp_count_matrix_pt->n_ij_e[pq]);
            dp_count_matrix_pt->n_ij_e[ij][bin_index] += dp_count_matrix_pt->n_ij_e[pq][bin_index_inner_loop];
          }
        }
      } /* end q-loop */
    }   /* end p-loop */
    new_c = MIN2(new_c, min_int);


    /* check for multibranch loops */
    //energy  = vrna_E_mb_loop_fast(fc, i, j, DMLi1, DMLi2);
    //new_c   = MIN2(new_c, energy);

    /* MULTI LOOP STRUCTURES */
    int min_ml = INF;
    if (!no_close) {
      /* dangle energies for multiloop closing stem */
      int tt    = rtype[type];
      int temp2 = fc->params->MLclosing;
      if (dangles == 2)
        temp2 += E_MLstem(tt, S1[j - 1], S1[i + 1], fc->params);
      else
        temp2 += E_MLstem(tt, -1, -1, fc->params);

      int u;
      for (u = i + TURN + 2; u < j - TURN - 2; u++) {
        int i1u   = fc->jindx[u] + (i+1);
        int u1j1  = fc->jindx[j-1] + (u+1);

        /* get distance to reference if closing the multiloop
         *  d = dbp(S_{i,j}, {i,j} + S_{i+1,u} + S_{u+1,j-1})
         */
        //if (!fc->matrices->fML[i1u])
        //  continue;

        //if (!fc->matrices->fM1[u1j1])
        //  continue;

        if ((fc->matrices->fML[i1u] != INF) && (fc->matrices->fM1[u1j1] != INF)) {
          int tmp_min_ml = MIN2(fc->matrices->c[ij],fc->matrices->fML[i1u]
                                                    + fc->matrices->fM1[u1j1]
                                                    + temp2
                                                    );
          min_ml = MIN2(min_ml, tmp_min_ml);

          fc->matrices->c[ij] = MIN2(fc->matrices->c[ij], fc->matrices->fML[i1u]  + fc->matrices->fM1[u1j1] + temp2);

          int joined_energy = fc->matrices->fML[i1u] + fc->matrices->fML[u1j1];
          int bin_index_m = get_energy_bin_index(min_energy, step_energy,energy_length, fc->matrices->fML[i1u]);
          int bin_index_m1 = get_energy_bin_index(min_energy, step_energy,energy_length, fc->matrices->fML[u1j1]);
          int bin_index_m_joined = get_energy_bin_index(min_energy, step_energy,energy_length, joined_energy);
          dp_count_matrix_pt->n_ij_e[ij][bin_index_m_joined] += dp_count_matrix_pt->n_ij_M_e[i1u][bin_index_m] * dp_count_matrix_pt->n_ij_M1_e[u1j1][bin_index_m1];   //dp_count_matrix_pt->n_ij_e[pq][bin_index_inner_loop];
        }
      }
    }

    new_c = MIN2(new_c, min_ml);


    if (dangle_model == 3) {
      /* coaxial stacking */
      energy  = vrna_E_mb_loop_stack(fc, i, j);
      new_c   = MIN2(new_c, energy);
    }

    /* remember stack energy for --noLP option */
    if (noLP) {
      stackEnergy = vrna_E_stack(fc, i, j);
      new_c       = MIN2(new_c, cc1[j - 1] + stackEnergy);
      cc[j]       = new_c;
      if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
        cc[j] -= fc->pscore[ij];

      e = cc1[j - 1] + stackEnergy;
    } else {
      e = new_c;
    }

    /* finally, check for auxiliary grammar rule(s) */
    if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux_c)) {
      energy  = fc->aux_grammar->cb_aux_c(fc, i, j, fc->aux_grammar->data);
      new_c   = MIN2(new_c, energy);
    }

    if (fc->type == VRNA_FC_TYPE_COMPARATIVE)
      e -= fc->pscore[ij];
  } /* end >> if (pair) << */

  return e;
}

int print_matrix(int *matrix, int *jiindx, int length_row, int length_col){
  int i,j;
  printf("\n");
  for(i=1; i <= length_row; i++){
    for(j=1; j <= length_col; j++){
      int mj = fmax(i,j);
      int mi = fmin(i,j);
      int ij = jiindx[mj] +mi;
      printf(" %7.4g", (float)matrix[ij]);
    }
    printf("\n");
  }

  return 0;
}

/* fill DP matrices */
PRIVATE int
fill_arrays(vrna_fold_compound_t *fc)
{
  int               i, j, ij, length, turn, uniq_ML, *indx, *f5, *c, *fML, *fM1;
  vrna_param_t      *P;
  vrna_mx_mfe_t     *matrices;
  vrna_ud_t         *domains_up;
  struct aux_arrays *helper_arrays;

  length      = (int)fc->length;
  indx        = fc->jindx;
  P           = fc->params;
  uniq_ML     = P->model_details.uniq_ML;
  turn        = P->model_details.min_loop_size;
  matrices    = fc->matrices;
  f5          = matrices->f5;
  c           = matrices->c;
  fML         = matrices->fML;
  fM1         = matrices->fM1;
  domains_up  = fc->domains_up;
  int dangles = P->model_details.dangles;
  short *S1            = fc->sequence_encoding;
  int circ = P->model_details.circ;


  int **n_ij_e = (int **)vrna_alloc(sizeof(int*) * ((length+1)*(length+1)));
  int **n_ij_A_e = (int **)vrna_alloc(sizeof(int*) * (length+1));
  //int **n_ij_B_e = (int **)vrna_alloc(sizeof(int*) * (length*length));
  int **n_ij_M_e = (int **)vrna_alloc(sizeof(int*) * ((length+1)*(length+1)));
  int **n_ij_M1_e = (int **)vrna_alloc(sizeof(int*) * ((length+1)*(length+1)));

  struct dp_counts_per_energy count_matrix_pt;
  count_matrix_pt.n_ij_e = n_ij_e;
  count_matrix_pt.n_ij_A_e = n_ij_A_e;
  //count_matrix_pt.n_ij_B_e = n_ij_B_e;
  count_matrix_pt.n_ij_M_e = n_ij_M_e;
  count_matrix_pt.n_ij_M1_e = n_ij_M1_e;


  int min_energy = -500; //(int)(vrna_mfe(fc,NULL) * 100); // TODO compute mfe.
  int max_energy = 1000; // TODO compute max_energy.
  int step_energy = 10; // 10 decakal. (smallest unit of energy computations)
  int range = max_energy - min_energy;
  int energy_length = (range/step_energy)+1;
  printf("min_energy: %d \n", min_energy);


  /* allocate memory for all helper arrays */
  helper_arrays = get_aux_arrays(length);

  /* pre-processing ligand binding production rule(s) */
  if (domains_up && domains_up->prod_cb)
    domains_up->prod_cb(fc, domains_up->data);

  /* prefill matrices with init contributions */
  for (j = 1; j <= length; j++)
    for (i = (j > turn ? (j - turn) : 1); i <= j; i++) {
      c[indx[j] + i] = fML[indx[j] + i] = INF;
      if (uniq_ML)
        fM1[indx[j] + i] = INF;
      printf(" %d", indx[j] + i);
    }
  printf("\n");

  /* start recursion */
  if (length <= turn) {
    /* clean up memory */
    free_aux_arrays(helper_arrays);

    /* return free energy of unfolded chain */
    return 0;
  }

  for (i = length - turn - 1; i >= 1; i--) {
    for (j = i + turn + 1; j <= length; j++) {
      ij = indx[j] + i;

      printf("alloc ij: %d %d %d\n", i,j, energy_length);

      int hc_decompose  = fc->hc->mx[length * i + j];
      int type  = fc->ptype[fc->jindx[j] + i];

      //prepare matrices (add third dimension)
      count_matrix_pt.n_ij_e[ij] = vrna_alloc(sizeof(int) * energy_length);
      //count_matrix_pt.n_ij_A_e[ij] = vrna_alloc(sizeof(int) * energy_length);
      //count_matrix_pt.n_ij_B_e[ij] = vrna_alloc(sizeof(int) * energy_length);
      count_matrix_pt.n_ij_M_e[ij] = vrna_alloc(sizeof(int) * energy_length);
      count_matrix_pt.n_ij_M1_e[ij] = vrna_alloc(sizeof(int) * energy_length);

      //printf("jiaddress from: %#08x \n", &(count_matrix_pt.n_ij_e[ij]));

      /* decompose subsegment [i, j] with pair (i, j) */
      //if (type)
      c[ij] = decompose_pair(fc, i, j, helper_arrays, min_energy, step_energy,energy_length, &count_matrix_pt);

      /* decompose subsegment [i, j] that is multibranch loop part with at least one branch */
      //fML[ij] = vrna_E_ml_stems_fast(fc, i, j, helper_arrays->Fmi, helper_arrays->DMLi);



      int temp2 = 0;
      printf("blub %d\n", i);
      /*
      if (dangles == 2)
        temp2 += E_MLstem(type, ((i > 1) || circ) ? S1[i - 1] : -1, ((j < length) || circ) ? S1[j + 1] : -1, P);
      else
        temp2 += E_MLstem(type, -1, -1, P);
      */
      if (dangles == 2)
        temp2 += E_MLstem(type, (i == 1) ? S1[length] : S1[i - 1], S1[j + 1], P);
      else
        temp2 += E_MLstem(type, -1, -1, P);


      /* now to the actual computations... */
      /* 1st E_M[ij] = E_M1[ij] = E_C[ij] + b */
      int energy_index;
     // if (matrices->c[ij]) {
            if (matrices->c[ij] != INF) {
              matrices->fML[ij] = matrices->fM1[ij] = temp2 + matrices->c[ij]; //TODO: compute correct energy bins.


              if(count_matrix_pt.n_ij_e[ij]){
              for(energy_index =0; energy_index < energy_length; energy_index++){

              count_matrix_pt.n_ij_M_e[ij][energy_index] = count_matrix_pt.n_ij_M1_e[ij][energy_index] = count_matrix_pt.n_ij_e[ij][energy_index];
//#ifdef COUNT_STATES
 //             matrices->N_M[ij][cnt1][cnt2 / 2] = matrices->N_M1[ij][cnt1][cnt2 / 2] = matrices->N_C[ij][cnt1][cnt2 / 2];
//#endif
            }
              }
        }
            else{
              matrices->fML[ij] = matrices->fM1[ij] = INF;
            }
      //}

      /* 2nd E_M[ij] = MIN(E_M[ij], E_M[i+1,j] + c) */
     // if (fc->matrices->fML[fc->jindx[j] + (i+1)]) {
            if (fc->matrices->fML[fc->jindx[j] + (i+1)] != INF) {
      matrices->fML[ij] = MIN2(matrices->fML[ij],matrices->fML[fc->jindx[j] + (i + 1)] + P->MLbase);
              //int energy_index;
              if(count_matrix_pt.n_ij_M_e[fc->jindx[j] + (i+1)]){
              for(energy_index =0; energy_index < energy_length; energy_index++){
                count_matrix_pt.n_ij_M_e[ij][energy_index] += count_matrix_pt.n_ij_M_e[fc->jindx[j] + (i+1)][energy_index];
//#ifdef COUNT_STATES
//                matrices->N_M[ij][cnt1 + dia][(cnt2 + dib) / 2] += matrices->N_M[my_iindx[i + 1] - j][cnt1][cnt2 / 2];
//#endif
              }
              }
            }
       //  }

      /* 3rd E_M[ij] = MIN(E_M[ij], E_M[i,j-1] + c) */
     // if (fc->matrices->fML[ij+1]) {
           if (fc->matrices->fML[fc->jindx[j-1] + i] != INF) {
      matrices->fML[ij] = MIN2(matrices->fML[ij],matrices->fML[fc->jindx[j-1] + i] + P->MLbase);
              //int energy_index;
              if(count_matrix_pt.n_ij_M_e[fc->jindx[j-1] + i]){
              for(energy_index =0; energy_index < energy_length; energy_index++){
                count_matrix_pt.n_ij_M_e[ij][energy_index] += count_matrix_pt.n_ij_M_e[fc->jindx[j-1] + i][energy_index];
//#ifdef COUNT_STATES
//                 matrices->N_M[ij][cnt1 + dja][(cnt2 + djb) / 2] += matrices->N_M[ij + 1][cnt1][cnt2 / 2];
//#endif
              }
              }
           }
       //   }

      /* 4th E_M1[ij] = MIN(E_M1[ij], E_M1[i,j-1] + c) */
     // if (fc->matrices->fM1[fc->jindx[j-1] +i]) {
            if (fc->matrices->fM1[fc->jindx[j-1] +i] != INF) {
     // if(count_matrix_pt.n_ij_M1_e[fc->jindx[j-1] +i]){
      matrices->fM1[ij] = MIN2(matrices->fM1[ij],matrices->fM1[fc->jindx[j-1] +i]+ P->MLbase);
              //int energy_index;
              if(count_matrix_pt.n_ij_M1_e[fc->jindx[j-1] +i]){
              for(energy_index =0; energy_index < energy_length; energy_index++){

                count_matrix_pt.n_ij_M1_e[ij][energy_index] += count_matrix_pt.n_ij_M1_e[fc->jindx[j-1] +i][energy_index];
//#ifdef COUNT_STATES
//                matrices->N_M1[ij][cnt1 + dja][(cnt2 + djb) / 2] += matrices->N_M1[ij + 1][cnt1][cnt2 / 2];
//#endif
              }
              }
      }
     //   }
      //}


      /* 5th E_M[ij] = MIN(E_M[ij], min(E_M[i,k] + E_M[k+1,j])) */
      if (j > TURN + 2) {
        int u;
        for (u = i + 1 + TURN; u <= j - 2 - TURN; u++) {
          //if (!fc->matrices->fML[fc->jindx[u]+i])
          //  continue;

          //if (!fc->matrices->fML[fc->jindx[j] + (u + 1)])
          //  continue;
          if ((fc->matrices->fML[fc->jindx[u]+i] != INF) && (fc->matrices->fML[fc->jindx[j] + (u + 1)] != INF))
            matrices->fML[ij] = MIN2(matrices->fML[ij], matrices->fML[fc->jindx[u]+i] + matrices->fML[fc->jindx[j] + (u + 1)]);
               //   if ((fc->matrices->fML[fc->jindx[u]+i] != INF) && (fc->matrices->fML[fc->jindx[j] + (u + 1)] != INF)) {
                    int energy_index;
                    if(count_matrix_pt.n_ij_M_e[fc->jindx[u]+i] != NULL && count_matrix_pt.n_ij_M1_e[fc->jindx[j] + (u + 1)] != NULL){
                    for(energy_index =0; energy_index < energy_length; energy_index++){
                      count_matrix_pt.n_ij_M_e[ij][energy_index] += count_matrix_pt.n_ij_M_e[fc->jindx[u]+i][energy_index] * count_matrix_pt.n_ij_M1_e[fc->jindx[j] + (u + 1)][energy_index];
                    }
                    }
//#ifdef COUNT_STATES
//                      matrices->N_M[ij][cnt1 + cnt3 + dia][(cnt2 + cnt4 + dib) / 2] += matrices->N_M[my_iindx[i] - u][cnt1][cnt2 / 2] * matrices->N_M1[my_iindx[u + 1] - j][cnt3][cnt4 / 2];
//#endif

                 // }

        }
      }


      /* decompose subsegment [i, j] that is multibranch loop part with exactly one branch */
      //if (uniq_ML)
      //  fM1[ij] = E_ml_rightmost_stem(i, j, fc);

      if ((fc->aux_grammar) && (fc->aux_grammar->cb_aux))
        fc->aux_grammar->cb_aux(fc, i, j, fc->aux_grammar->data);
    } /* end of j-loop */

    rotate_aux_arrays(helper_arrays, length);

    printf("\nPairMatrix:");
    print_matrix(c, fc->jindx, length, length);
    printf("\nM_Matrix:");
    print_matrix(fML, fc->jindx, length, length);
    printf("\nM1_Matrix:");
    print_matrix(fM1, fc->jindx, length, length);
    printf("\nF_Matrix:");
  //  print_matrix(f5, fc->jindx, 1, length+1);

  } /* end of i-loop */
    /* calculate energies of 5' fragments */
 // (void)vrna_E_ext_loop_5(fc);

  int x;
  for(x = 0; x <= length; x++){
    count_matrix_pt.n_ij_A_e[x] = vrna_alloc(sizeof(int) * energy_length);
  }


  for (j = TURN + 2; j <= length; j++) {
    int energy_index;
    /* j-1 is unpaired ... */
    if(count_matrix_pt.n_ij_A_e[j]){
    for(energy_index =0; energy_index < energy_length; energy_index++){
      count_matrix_pt.n_ij_A_e[j][energy_index] += count_matrix_pt.n_ij_A_e[j-1][energy_index];
    }
    }
   // matrices->N_F5[j][cnt1 + da][(cnt2 + db) / 2] += matrices->N_F5[j - 1][cnt1][cnt2 / 2];

    /* j pairs with 1 */
    if(count_matrix_pt.n_ij_e[fc->jindx[j] + 1]){
    for(energy_index =0; energy_index < energy_length; energy_index++){
          count_matrix_pt.n_ij_A_e[j][energy_index] += count_matrix_pt.n_ij_e[fc->jindx[j] + 1][energy_index];
        }
       // matrices->N_F5[j][cnt1][cnt2 / 2] += matrices->N_C[my_iindx[1] - j][cnt1][cnt2 / 2];
    }


    /* j pairs with some other nucleotide -> see below */
        for (i = j - TURN - 1; i > 1; i--) {
          ij    = fc->jindx[j] + i;
          int type  = fc->ptype[ij];
          if (type) {

          //  if (!fc->matrices->c[ij])
          //    continue;

                    if (fc->matrices->f5[i - 1] != INF && fc->matrices->c[ij] != INF) {
                      for(energy_index =0; energy_index < energy_length; energy_index++){
                        count_matrix_pt.n_ij_A_e[j][energy_index] += count_matrix_pt.n_ij_A_e[i-1][energy_index] * count_matrix_pt.n_ij_e[ij][energy_index];
   // #ifdef COUNT_STATES
    //                    matrices->N_F5[j][cnt1 + cnt3 + d1a][(cnt2 + cnt4 + d1b) / 2] += matrices->N_F5[i - 1][cnt3][cnt4 / 2] * matrices->N_C[ij][cnt1][cnt2 / 2];
    //#endif
                      }

          }
        }
        }
  }

  printf("\n");
  int fidx = 0;
  for(; fidx < length+1; fidx++){
    printf(" %d", matrices->f5[fidx]);
  }
  printf("\n");

 // print_matrix(f5, fc->jindx, 1, length+1);


  int energy_index = 0;
  for(energy_index = 0; energy_index < energy_length; energy_index++){
    if(count_matrix_pt.n_ij_A_e[length]){
      int counts = count_matrix_pt.n_ij_A_e[length][energy_index];
      int energy = min_energy + energy_index * step_energy;
      printf("energy: %d, counts: %d \n", energy, counts);
    }
  }


  /* clean up memory */
  free_aux_arrays(helper_arrays);

  for (i = length - turn - 1; i >= 1; i--) {
      for (j = i + turn + 1; j <= length; j++) {
        ij = indx[j] + i;
        free(count_matrix_pt.n_ij_e[ij]);
        //free(count_matrix_pt.n_ij_A_e[ij]);
        //free(count_matrix_pt.n_ij_B_e[ij]);
        free(count_matrix_pt.n_ij_M_e[ij]);
        free(count_matrix_pt.n_ij_M1_e[ij]);
      }
  }

  for(i = 0; i <= length; i++){
    free(count_matrix_pt.n_ij_A_e[i]);
  }

  free(count_matrix_pt.n_ij_e);
  free(count_matrix_pt.n_ij_A_e);
  //free(count_matrix_pt.n_ij_B_e);
  free(count_matrix_pt.n_ij_M_e);
  free(count_matrix_pt.n_ij_M1_e);


  return f5[length];
}


int main(){
printf("hi\n");
char *sequence = "ACGUACGUUGCA";
vrna_fold_compound_t  *fc;
vrna_md_t md;
set_model_details(&md);
md.uniq_ML =1;
md.noLP = 0;
md.circ = 0;
md.dangles =2;
//md.temperature  = 37;
vrna_param_t *params               = vrna_params(&md);
//params->model_details.circ = is_circular;
fc = vrna_fold_compound(sequence, &(params->model_details), VRNA_OPTION_DEFAULT | VRNA_OPTION_MFE);

if (!vrna_fold_compound_prepare(fc, VRNA_OPTION_MFE)) {
  vrna_message_warning("vrna_mfe@mfe.c: Failed to prepare vrna_fold_compound");
}

int count = fill_arrays(fc);

/*
    vrna_mfe(fc,NULL);
int length = fc->length;
printf("\nPairMatrix:");
print_matrix(fc->matrices->c, fc->jindx, length, length);
printf("\nM_Matrix:");
print_matrix(fc->matrices->fML, fc->jindx, length, length);
printf("\nM1_Matrix:");
print_matrix(fc->matrices->fM1, fc->jindx, length, length);
printf("\nF_Matrix:");
print_matrix(fc->matrices->f5, fc->jindx, 1, length+1);
*/
printf("count: %d \n", count);

free(params);
vrna_fold_compound_free(fc);

}
