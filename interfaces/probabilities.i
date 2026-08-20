/**********************************************/
/* BEGIN interface for Equilibrium            */
/* Probabilities, i.e. outside Partition      */
/* Function computation                       */
/**********************************************/

%ignore export_bppm;
%ignore mean_bp_distance_pr;
%ignore stackProb;
%ignore mean_bp_dist;
%ignore assign_plist_gquad_from_pr;

%extend vrna_fold_compound_t{

  double
  mean_bp_distance()
  {
    return vrna_mean_bp_distance($self);
  }

  double
  ensemble_defect(std::string structure,
                  unsigned int options = VRNA_BRACKETS_RND)
  {
    double ed;
    short int         *pt;

    pt = vrna_ptable_from_string(structure.c_str(), options);

    ed = vrna_ensemble_defect_pt($self, pt);

    free(pt);

    return ed;
  }

  double
  ensemble_defect(std::vector<int> pair_table)
  {
    std::vector<short> pt_v_short;
    transform(pair_table.begin(), pair_table.end(), back_inserter(pt_v_short), convert_vecint2vecshort);
    return vrna_ensemble_defect_pt($self, (short*)&pt_v_short[0]);
  }

  double
  ensemble_defect(var_array<short> const &pt)
  {
    return vrna_ensemble_defect_pt($self, pt.data);
  }

  std::vector<double>
  positional_entropy(void)
  {
    unsigned int        n;
    double              *pos_ent;
    std::vector<double> dv;

    n       = $self->length;
    pos_ent = vrna_positional_entropy($self);

    if (pos_ent)
      dv.assign(pos_ent, pos_ent + (n + 1));

    return dv;
  }

  double
  pr_structure(std::string structure)
  {
    return vrna_pr_structure($self, structure.c_str());
  }

  double
  pr_energy(double e)
  {
    return vrna_pr_energy($self, e);
  }


  std::vector<std::vector<double> >
  bpp(bool symmetric = false)
  {
    std::vector<std::vector<double> > probabilities;
    vrna_fold_compound_t *vc = $self;
    if(vc->exp_matrices && vc->exp_matrices->probs){
      int turn, i, j, *idx, n;
      FLT_OR_DBL *probs;

      n     = vc->length;
      idx   = vc->iindx;
      probs = vc->exp_matrices->probs;

      for(i = 0; i <= n; i++)
        probabilities.push_back(std::vector<double>(n + 1, 0.));

      for(i = 1; i <= n; i++) {
        for(j = i + 1; j <= n; j++) {
          probabilities[i][j] = (double)probs[idx[i] - j];
          if (symmetric)
            probabilities[j][i] = (double)probs[idx[i] - j];
        }
      }
    }

    return probabilities;
  }


  std::vector<vrna_ep_t>
  stack_prob(double cutoff = 1e-5)
  {
    std::vector<vrna_ep_t > ep_v;
    vrna_ep_t               *ptr, *plist;

    plist = vrna_stack_prob($self, cutoff);

    for (ptr = plist; ptr->i && ptr->j; ptr++) {
      vrna_ep_t pl;
      pl.i = ptr->i;
      pl.j = ptr->j;
      pl.p = ptr->p;
      pl.type = ptr->type;
      ep_v.push_back(pl);
    }

    free(plist);

    return ep_v;
  }
}


/* scripting language access through 'fold_compound' instead of 'vrna_fold_compound_t' */
%rename(eqm_probs_opt) vrna_eqm_probs_opt_t;

#ifdef SWIGPYTHON
/* create object oriented interface for vrna_eqm_probs_opt_t */
%extend vrna_eqm_probs_opt_t {

 std::string
  __str__()
  {
    std::ostringstream out;

    out << "{ mode: " ;
    if ($self->mode & VRNA_EQM_PROBS_BPP)
      out << "EQM_PROBS_BPP";
    out << ", up_length: " << $self->up_length;
    out << ", hx_length: " << $self->hx_length;
    out << " }";

    return std::string(out.str());
  }

%pythoncode %{
def __repr__(self):
    # reformat string representation (self.__str__()) to something
    # that looks like a constructor argument list
    strthis = self.__str__().replace(": ", "=").replace("{ ", "").replace(" }", "")
    return  "%s.%s(%s)" % (self.__class__.__module__, self.__class__.__name__, strthis) 
%}
}
#endif


%rename (eqm_probs_up)  my_eqm_probs_up;

%{
    std::vector<vrna_ep_t>
    my_eqm_probs_up(vrna_fold_compound_t *fc,
                    unsigned int         min_u,
                    unsigned int         max_u,
                    double               cutoff)
    {
      std::vector<vrna_ep_t> result;

      vrna_array(vrna_ep_t) c_result;
      c_result = vrna_eqm_probs_up(fc, min_u, max_u, cutoff);

      result = {&(c_result[0]), &(c_result[0]) + vrna_array_size(c_result) };

      vrna_array_free(c_result);

      return result;
    }
%}

%extend vrna_fold_compound_t {

    std::vector<vrna_ep_t>
    my_eqm_probs_up(unsigned int         min_u,
                    unsigned int         max_u,
                    double               cutoff)
    {
      std::vector<vrna_ep_t> result;

      vrna_array(vrna_ep_t) c_result;
      c_result = vrna_eqm_probs_up($self, min_u, max_u, cutoff);

      result = {&(c_result[0]), &(c_result[0]) + vrna_array_size(c_result) };

      vrna_array_free(c_result);

      return result;
    }
}

std::vector<vrna_ep_t>
my_eqm_probs_up(vrna_fold_compound_t *fc,
                unsigned int         min_u,
                unsigned int         max_u,
                double               cutoff);


%constant unsigned int LOOP_TYPE_EXT  = VRNA_LOOP_TYPE_EXT;
%constant unsigned int LOOP_TYPE_HP   = VRNA_LOOP_TYPE_HP;
%constant unsigned int LOOP_TYPE_INT  = VRNA_LOOP_TYPE_INT;
%constant unsigned int LOOP_TYPE_MB   = VRNA_LOOP_TYPE_MB;
%constant unsigned int LOOP_TYPE_ANY  = VRNA_LOOP_TYPE_ANY;

%constant unsigned int EQM_PROBS_BPP            = VRNA_EQM_PROBS_BPP;
%constant unsigned int EQM_PROBS_UP             = VRNA_EQM_PROBS_UP;
%constant unsigned int EQM_PROBS_HX             = VRNA_EQM_PROBS_HX;
%constant unsigned int EQM_PROBS_SPLIT          = VRNA_EQM_PROBS_SPLIT;
%constant unsigned int EQM_PROBS_STORE_OUTSIDE  = VRNA_EQM_PROBS_STORE_OUTSIDE;
%constant unsigned int EQM_PROBS_OUTSIDE_PF     = VRNA_EQM_PROBS_OUTSIDE_PF;


%include  <ViennaRNA/probabilities/basic.h>
%include  <ViennaRNA/probabilities/basepairs.h>
%include  <ViennaRNA/probabilities/structures.h>

%{
double
get_pr(int i,
       int j)
{
  int ii;
  if (i>j) {ii=i; i=j; j=ii;}
  return pr[iindx[i]-j];
}
%}
double get_pr(int i, int j);
/* Get probability of pair i.j from the pr array */
