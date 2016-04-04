/**********************************************/
/* BEGIN interface for energy evaluation      */
/**********************************************/


%extend vrna_fold_compound_t{

 float eval_structure(const char *structure){
	  return vrna_eval_structure($self,structure);
  }
  /*calculate MFE of given pairtable*/
  float eval_structure_pt(std::vector<int> pt)
  {
	  std::vector<short> vc;
	  transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);
	  return vrna_eval_structure_pt($self,(const short*)&vc[0]);
  }
  
  
  /*MFE of given structure, but now with different FileHandler for verbose, NULL = STDOUT*/
  float eval_structure_verbose(char *structure, FILE *file)
  {
	  return vrna_eval_structure_verbose($self,structure,file);
  }
  
 /*MFE of given pairtable, with different FileHandler for verbose, Default value = NULL + STDOUT*/
  float eval_structure_pt_verbose(std::vector<int> pt, FILE *file)
  {
	  std::vector<short> vc;
	  transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);
	  return vrna_eval_structure_pt_verbose($self,(const short*)&vc[0],file);
  }
  
  /*returns the energy for a structure to a given set of alignment sequences*/
  float eval_covar_structure2(char * structure)
  {    
	  return vrna_eval_covar_structure($self, structure);
  }
  

  
  /*returns the energy of a loop specified by i to pt[i]*/
  float eval_loop_pt(int i, std::vector<int> pt)
  {
	  std::vector<short> vc;
	  transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);
	  return vrna_eval_loop_pt($self,i,(const short*)&vc[0]);
  }
  
  /*returns the energy change by introducing a move on a given structure*/
  float eval_move(const char *structure,int m1, int m2)
  {
	  return vrna_eval_move($self,structure,m1,m2);
  }
  /*returns the energy change by introducing a move on a given pairtable*/
  float eval_move_pt(std::vector<int> pt,int m1, int m2)
  {
	  std::vector<short> vc;
	  transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);
	  return vrna_eval_move_pt($self,((short*)&vc[0]),m1,m2);   /*attention here no cosnt short* as argument*/
  }
}


%ignore energy_of_struct_par;
%ignore energy_of_circ_struct_par;
%ignore energy_of_gquad_struct_par;
%ignore energy_of_struct_pt_par;

%include  <ViennaRNA/eval.h>
