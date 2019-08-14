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
    return vrna_eval_structure_pt($self, (const short*)&vc[0]);
  }
  
  /*  compute free energy for structure given in dot-bracket notation, 
      but now with different FileHandler for verbose, NULL = STDOUT */
  float eval_structure_verbose(char *structure, FILE *nullfile = NULL){

    return vrna_eval_structure_verbose($self, structure, nullfile);
  }
  
  /* compute free energy for structure given in pairtable (verbose), with different FileHandler for verbose, Default value = NULL + STDOUT*/
  int eval_structure_pt_verbose(std::vector<int> pt, FILE *nullfile = NULL){

    std::vector<short> vc;
    transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);
    return vrna_eval_structure_pt_verbose($self, (const short*)&vc[0], nullfile);
  }
  
  /* compute covariance contributions for consensus structure given in dot-bracket notation */
  float eval_covar_structure(char * structure){

    return vrna_eval_covar_structure($self, structure);
  }

  /* returns the energy of a loop specified by i to pt[i] */
  int eval_loop_pt(int i, std::vector<int> pt){

    std::vector<short> vc;
    transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);
    return vrna_eval_loop_pt($self, i, (const short*)&vc[0]);
  }

  /* returns the energy change by introducing a move on a given structure */
  float eval_move(const char *structure, int m1, int m2){

    return vrna_eval_move($self, structure, m1, m2);
  }

  /* returns the energy change by introducing a move on a given pairtable */
  int eval_move_pt(std::vector<int> pt, int m1, int m2){

    std::vector<short> vc;
    transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);
    return vrna_eval_move_pt($self, ((short*)&vc[0]), m1, m2);   /*attention here no const short* as argument*/
  }
}

%ignore energy_of_struct_par;
%ignore energy_of_circ_struct_par;
%ignore energy_of_gquad_struct_par;
%ignore energy_of_struct_pt_par;


/*
  Below follow the wrappers for simplified energy evaluations
  Here, we create a single overloaded wrapper that allows to:
  1. use single sequences or sequence alignments as first parameter
  2. omit file handle (last parameter)
  3. omit verbosity level (third parameter)

  Thus, we do not wrap all variations of a particular function
  but provide a single one to rule them all!
 */
%rename (eval_structure_simple) my_eval_structure_simple;
%rename (eval_circ_structure) my_eval_circ_structure;
%rename (eval_gquad_structure) my_eval_gquad_structure;
%rename (eval_circ_gquad_structure) my_eval_circ_gquad_structure;
%rename (eval_structure_pt_simple) my_eval_structure_pt_simple;


%{

  float
  my_eval_structure_simple(std::string sequence,
                           std::string structure,
                           int         verbosity_level = VRNA_VERBOSITY_QUIET,
                           FILE        *file = NULL)
  {
    return vrna_eval_structure_simple_v(sequence.c_str(), structure.c_str(), verbosity_level, file);
  }

  float
  my_eval_circ_structure(std::string sequence,
                         std::string structure,
                         int         verbosity_level = VRNA_VERBOSITY_QUIET,
                         FILE        *file = NULL)
  {
    return vrna_eval_circ_structure_v(sequence.c_str(), structure.c_str(), verbosity_level, file);
  }

  float
  my_eval_gquad_structure(std::string sequence,
                          std::string structure,
                          int         verbosity_level = VRNA_VERBOSITY_QUIET,
                          FILE        *file = NULL)
  {
    return vrna_eval_gquad_structure_v(sequence.c_str(), structure.c_str(), verbosity_level, file);
  }

  float
  my_eval_circ_gquad_structure(std::string sequence,
                               std::string structure,
                               int         verbosity_level = VRNA_VERBOSITY_QUIET,
                               FILE        *file = NULL)
  {
    return vrna_eval_circ_gquad_structure_v(sequence.c_str(), structure.c_str(), verbosity_level, file);
  }

  float
  my_eval_structure_simple(std::vector<std::string> alignment,
                           std::string structure,
                           int         verbosity_level = VRNA_VERBOSITY_QUIET,
                           FILE        *file = NULL)
  {
    std::vector<const char*>  vc;

    std::transform(alignment.begin(), alignment.end(), std::back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */

    return vrna_eval_consensus_structure_simple_v((const char **)&vc[0], structure.c_str(), verbosity_level, file);
  }

  float
  my_eval_circ_structure(std::vector<std::string> alignment,
                         std::string structure,
                         int         verbosity_level = VRNA_VERBOSITY_QUIET,
                         FILE        *file = NULL)
  {
    std::vector<const char*>  vc;

    std::transform(alignment.begin(), alignment.end(), std::back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */

    return vrna_eval_circ_consensus_structure_v((const char **)&vc[0], structure.c_str(), verbosity_level, file);
  }

  float
  my_eval_gquad_structure(std::vector<std::string> alignment,
                          std::string structure,
                          int         verbosity_level = VRNA_VERBOSITY_QUIET,
                          FILE        *file = NULL)
  {
    std::vector<const char*>  vc;

    std::transform(alignment.begin(), alignment.end(), std::back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */

    return vrna_eval_gquad_consensus_structure_v((const char **)&vc[0], structure.c_str(), verbosity_level, file);
  }

  float
  my_eval_circ_gquad_structure(std::vector<std::string> alignment,
                               std::string structure,
                               int         verbosity_level = VRNA_VERBOSITY_QUIET,
                               FILE        *file = NULL)
  {
    std::vector<const char*>  vc;

    std::transform(alignment.begin(), alignment.end(), std::back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */

    return vrna_eval_circ_gquad_consensus_structure_v((const char **)&vc[0], structure.c_str(), verbosity_level, file);
  }

  float
  my_eval_structure_pt_simple(std::string sequence,
                              std::vector<int> pt,
                              int         verbosity_level = VRNA_VERBOSITY_QUIET,
                              FILE        *file = NULL)
  {
    std::vector<short> vc;

    transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);

    return vrna_eval_structure_pt_simple_v(sequence.c_str(), (const short*)&vc[0], verbosity_level, file);
  }

  float
  my_eval_structure_pt_simple(std::vector<std::string> alignment,
                              std::vector<int> pt,
                              int         verbosity_level = VRNA_VERBOSITY_QUIET,
                              FILE        *file = NULL)
  {
    std::vector<const char*>  vc;
    std::vector<short> ptv;

    std::transform(alignment.begin(), alignment.end(), std::back_inserter(vc), convert_vecstring2veccharcp);
    vc.push_back(NULL); /* mark end of sequences */

    transform(pt.begin(), pt.end(), back_inserter(ptv), convert_vecint2vecshort);

    return vrna_eval_consensus_structure_pt_simple_v((const char **)&vc[0], (const short*)&ptv[0], verbosity_level, file);
  }

%}

float
my_eval_structure_simple(std::string sequence,
                         std::string structure,
                           int         verbosity_level = VRNA_VERBOSITY_QUIET,
                           FILE        *file = NULL);

float
my_eval_circ_structure(std::string sequence,
                       std::string structure,
                       int         verbosity_level = VRNA_VERBOSITY_QUIET,
                       FILE        *file = NULL);

float
my_eval_gquad_structure(std::string sequence,
                        std::string structure,
                        int         verbosity_level = VRNA_VERBOSITY_QUIET,
                        FILE        *file = NULL);

float
my_eval_circ_gquad_structure(std::string sequence,
                             std::string structure,
                             int         verbosity_level = VRNA_VERBOSITY_QUIET,
                             FILE        *file = NULL);

float
my_eval_structure_simple(std::vector<std::string> alignment,
                         std::string structure,
                           int         verbosity_level = VRNA_VERBOSITY_QUIET,
                           FILE        *file = NULL);

float
my_eval_circ_structure(std::vector<std::string> alignment,
                       std::string structure,
                       int         verbosity_level = VRNA_VERBOSITY_QUIET,
                       FILE        *file = NULL);

float
my_eval_gquad_structure(std::vector<std::string> alignment,
                        std::string structure,
                        int         verbosity_level = VRNA_VERBOSITY_QUIET,
                        FILE        *file = NULL);

float
my_eval_circ_gquad_structure(std::vector<std::string> alignment,
                             std::string structure,
                             int         verbosity_level = VRNA_VERBOSITY_QUIET,
                             FILE        *file = NULL);

float
my_eval_structure_pt_simple(std::string sequence,
                            std::vector<int> pt,
                            int         verbosity_level = VRNA_VERBOSITY_QUIET,
                            FILE        *file = NULL);


float
my_eval_structure_pt_simple(std::vector<std::string> alignment,
                            std::vector<int> pt,
                            int         verbosity_level = VRNA_VERBOSITY_QUIET,
                            FILE        *file = NULL);


%ignore vrna_eval_structure_simple;
%ignore vrna_eval_structure_simple_v;
%ignore vrna_eval_structure_simple_verbose;
%ignore vrna_eval_circ_structure;
%ignore vrna_eval_circ_structure_v;
%ignore vrna_eval_gquad_structure;
%ignore vrna_eval_gquad_structure_v;
%ignore vrna_eval_circ_gquad_structure;
%ignore vrna_eval_circ_gquad_structure_v;
%ignore vrna_eval_consensus_structure_simple;
%ignore vrna_eval_consensus_structure_simple_verbose;
%ignore vrna_eval_consensus_structure_simple_v;
%ignore vrna_eval_circ_consensus_structure;
%ignore vrna_eval_circ_consensus_structure_v;
%ignore vrna_eval_gquad_consensus_structure;
%ignore vrna_eval_gquad_consensus_structure_v;
%ignore vrna_eval_circ_gquad_consensus_structure;
%ignore vrna_eval_circ_gquad_consensus_structure_v;
%ignore vrna_eval_structure_pt_simple;
%ignore vrna_eval_consensus_structure_pt_simple;


%include  <ViennaRNA/eval.h>
