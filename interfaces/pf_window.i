/**********************************************/
/* BEGIN interface for sliding-window         */
/* partition function                         */
/**********************************************/

%ignore update_pf_paramsLP;
%ignore update_pf_paramsLP_par;
%ignore pfl_fold;
%ignore pfl_fold_par;
%ignore putoutpU_prob_par;
%ignore putoutpU_prob;
%ignore putoutpU_prob_bin_par;
%ignore putoutpU_prob_bin;
%ignore init_pf_foldLP;

%rename (pfl_fold) my_pfl_fold;

%{
#include <vector>

  std::vector<vrna_ep_t> my_pfl_fold(std::string sequence, int w, int L, double cutoff)
  {
    std::vector<vrna_ep_t > vplist;
    vrna_ep_t *ptr, *plist;

    plist = vrna_pfl_fold(sequence.c_str(), w, L, (float)cutoff);

    for (ptr = plist; ptr->i && ptr->j; ptr++) {
      vrna_ep_t pl;
      pl.i    = ptr->i;
      pl.j    = ptr->j;
      pl.p    = ptr->p;
      pl.type = ptr->type;
      vplist.push_back(pl);
    }
    free(plist);

    return vplist;
  }

  std::vector<std::vector<double> > pfl_fold_up(std::string sequence,
                                                int ulength,
                                                int window_size,
                                                int max_bp_span)
  {
    double **up = vrna_pfl_fold_up(sequence.c_str(), ulength, window_size, max_bp_span);

    std::vector<std::vector<double> > up_vec;

    std::vector<double> nullvec(ulength + 1, 0.);

    /* insert a 0th element, since we start a 1-based N x M matrix here */
    up_vec.push_back(nullvec);
    free(up[0]);
    for (unsigned int i = 1; i <= sequence.length(); i++) {
      std::vector<double> row;
      /* insert a 0th element, again, everything should be 1-based */
      row.push_back(0.);

      /* add remaining elements for this row */
      for (int j = 1; j <= ulength; j++) {
        row.push_back(up[i][j]);
      }

      /* free memory of i-th row in up array */
      free(up[i]);

      up_vec.push_back(row);
    }
    free(up);

    return up_vec;
  }
%}

#ifdef SWIGPYTHON
%feature("autodoc") my_pfl_fold;
%feature("kwargs") my_pfl_fold;
%feature("autodoc") pfl_fold_up;
%feature("kwargs") pfl_fold_up;
#endif

std::vector<vrna_ep_t> my_pfl_fold(std::string sequence, int w, int L, double cutoff);
std::vector<std::vector<double> > pfl_fold_up(std::string sequence,
                                              int ulength,
                                              int window_size,
                                              int max_bp_span);

%constant unsigned int EXT_LOOP = VRNA_EXT_LOOP;
%constant unsigned int HP_LOOP  = VRNA_HP_LOOP;
%constant unsigned int INT_LOOP = VRNA_INT_LOOP;
%constant unsigned int MB_LOOP  = VRNA_MB_LOOP;
%constant unsigned int ANY_LOOP = VRNA_ANY_LOOP;

%constant unsigned int PROBS_WINDOW_BPP       = VRNA_PROBS_WINDOW_BPP;
%constant unsigned int PROBS_WINDOW_UP        = VRNA_PROBS_WINDOW_UP;
%constant unsigned int PROBS_WINDOW_STACKP    = VRNA_PROBS_WINDOW_STACKP;
%constant unsigned int PROBS_WINDOW_UP_SPLIT  = VRNA_PROBS_WINDOW_UP_SPLIT;
%constant unsigned int PROBS_WINDOW_PF        = VRNA_PROBS_WINDOW_PF;

%include  <ViennaRNA/LPfold.h>
