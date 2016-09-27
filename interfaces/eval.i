/**********************************************/
/* BEGIN interface for energy evaluation      */
/**********************************************/


%extend vrna_fold_compound_t{

  float eval_structure(const char *structure){

    return vrna_eval_structure($self,structure);
  }

  /* calculate free energy for structure given in pairtable */
  int eval_structure_pt(std::vector<int> pt){

    std::vector<short> vc;
    transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);
    return vrna_eval_structure_pt($self,(const short*)&vc[0]);
  }
  
  /*  compute free energy for structure given in dot-bracket notation, 
      but now with different FileHandler for verbose, NULL = STDOUT */
  float eval_structure_verbose(char *structure, FILE *file){

    return vrna_eval_structure_verbose($self,structure,file);
  }
  
  /* compute free energy for structure given in pairtable (verbose), with different FileHandler for verbose, Default value = NULL + STDOUT*/
  int eval_structure_pt_verbose(std::vector<int> pt, FILE *file){

    std::vector<short> vc;
    transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);
    return vrna_eval_structure_pt_verbose($self,(const short*)&vc[0],file);
  }
  
  /* compute covariance contributions for consensus structure given in dot-bracket notation */
  float eval_covar_structure(char * structure){

    return vrna_eval_covar_structure($self, structure);
  }

  /* returns the energy of a loop specified by i to pt[i] */
  int eval_loop_pt(int i, std::vector<int> pt){

    std::vector<short> vc;
    transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);
    return vrna_eval_loop_pt($self,i,(const short*)&vc[0]);
  }

  /* returns the energy change by introducing a move on a given structure */
  float eval_move(const char *structure,int m1, int m2){

    return vrna_eval_move($self,structure,m1,m2);
  }

  /* returns the energy change by introducing a move on a given pairtable */
  int eval_move_pt(std::vector<int> pt,int m1, int m2){

    std::vector<short> vc;
    transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);
    return vrna_eval_move_pt($self,((short*)&vc[0]),m1,m2);   /*attention here no const short* as argument*/
  }
}

%ignore energy_of_struct_par;
%ignore energy_of_circ_struct_par;
%ignore energy_of_gquad_struct_par;
%ignore energy_of_struct_pt_par;

%include  <ViennaRNA/eval.h>


%extend vrna_fold_compound_t{

  int E_ext_loop(int i, int j){
    return E_ext_loop(i, j, $self);
  }

  int eval_hp_loop(int i, int j){
    return vrna_eval_hp_loop($self, i, j);
  }

  int eval_int_loop(int i, int j, int k, int l){
    return vrna_eval_int_loop($self, i, j, k, l);
  }

}

%include  <ViennaRNA/exterior_loops.h>
%include  <ViennaRNA/hairpin_loops.h>
%include  <ViennaRNA/interior_loops.h>
%include  <ViennaRNA/multibranch_loops.h>
