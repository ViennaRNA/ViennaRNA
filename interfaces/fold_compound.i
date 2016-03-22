/**********************************************/
/* BEGIN interface for fold compound          */
/**********************************************/

/* ignore all data structures we handle elsewhere */
%ignore PAIR;
%ignore plist;
%ignore cpair;
%ignore sect;
%ignore bondT;
%ignore pu_contrib;
%ignore interact;
%ignore pu_out;
%ignore constrain;
//%ignore duplexT;
%ignore folden;
%ignore snoopT;
%ignore dupVar;

/* add callback binding methods for fold_compound */
%include callbacks-fc.i
%include callbacks-sc.i


/* start constructing a sane interface to vrna_fold_compound_t */

%rename(fc_type) vrna_fc_type_e;

/* scripting language access through 'fold_compound' instead of 'vrna_fold_compound_t' */
%rename(fold_compound) vrna_fold_compound_t;

/* no default constructor / destructor */
%nodefaultctor vrna_fold_compound_t;
%nodefaultdtor vrna_fold_compound_t;


/* hide all attributes in vrna_fold_compound_t */
typedef struct {} vrna_fold_compound_t;

/* scripting language takes ownership of objects returned by mfe() method */
%newobject vrna_fold_compound_t::mfe;

/* create object oriented interface for vrna_fold_compount_t */
%extend vrna_fold_compound_t {

  
  
  /* the default constructor, *md is optional, for single sequences*/
  vrna_fold_compound_t(const char *sequence, vrna_md_t *md=NULL, unsigned int options=VRNA_OPTION_MFE){
    return vrna_fold_compound(sequence, md, options);
  }
   /*the constructor for alignments, *md is optional  
   vrna_fold_compound_t(const char **sequences, vrna_md_t *md=NULL, unsigned int options=VRNA_OPTION_MFE){
   return vrna_fold_compound_comparative(sequences, md, options);
  }*/
  /*vrna_fold_compound_t(const char **sequences){
    return vrna_fold_compound_comparative(sequences, NULL, VRNA_OPTION_MFE);
  }*/
  /* constructor for distance class partitioning, *md is optional, for single sequences
  vrna_fold_compound_t(const char *sequence,char *s1,char *s2, vrna_md_t *md=NULL, unsigned int options=VRNA_OPTION_MFE){
    return vrna_fold_compound_TwoD(sequence,s1,s2, md, options);
  }*/
  
  
  
  
  ~vrna_fold_compound_t(){
    vrna_fold_compound_free($self);
  }

  vrna_fc_type_e type(){
    return $self->type;
  }

  char *mfe(float *OUTPUT){
    char *structure = (char *)vrna_alloc(sizeof(char) * ($self->length + 1));
    *OUTPUT = vrna_mfe($self, structure);
    return structure;
  }
  /*MFE for 2 dimerisated RNA strands*/
  char *mfe_dimer(float *OUTPUT){
    char *structure = vrna_alloc(sizeof(char) * ($self->length + 1));
    *OUTPUT = vrna_mfe_dimer($self, structure);
    return structure;
  }
  
  char* testFunction(const char** str,vrna_md_t *md=NULL,unsigned int options=VRNA_OPTION_MFE)
  {
	  return "geht eh";
  }
  char* testFunction2(const char** str,unsigned int options=VRNA_OPTION_MFE,vrna_md_t *md=NULL)
  {
	  return "geht auch";
  }
  
  /*calculate MFE of given structure in dot bracket format*/
  
  
  /*####################
   * from eval.h
#####################################*/
  
  
  float eval_structure(const char *structure){
	  return vrna_eval_structure($self,structure);
  }
  /*calculate MFE of given pairtable*/
  float eval_structure_pt(short * pt)
  {	  
	return vrna_eval_structure_pt($self,pt);
	  
}
  
  /*MFE of given structure, but now with different FileHandler for verbose, Default value = NULL + STDOUT*/
  float eval_structure_verbose(const char *structure, FILE *file=NULL)
  {
	  return vrna_eval_structure_verbose($self,structure,file);
  }
 /*MFE of given pairtable, with different FileHandler for verbose, Default value = NULL + STDOUT*/
  float eval_structure_pt_verbose(short *pt, FILE *file=NULL)
  {
	  return vrna_eval_structure_pt_verbose($self,pt,file);
  }
  
  /*returns the energy for a secondary consensus structure calculated for a given set of alignment sequences*/
  float eval_covar_structure(char *structure)
  {
	  return vrna_eval_covar_structure($self,structure);
  }

  
  /*returns the energy of a loop specified by i to pt[i]*/
  float eval_loop_pt(int i, short *pt)
  {
	  return vrna_eval_loop_pt($self,i,pt);
  }
  
  /*returns the energy change by introducing a move on a given structure*/
  float eval_move(const char *structure,int m1, int m2)
  {
	  return vrna_eval_move($self,structure,m1,m2);
  }
  /*returns the energy change by introducing a move on a given pairtable*/
  float eval_move_pt(short *pt,int m1, int m2)
  {
	  return vrna_eval_move_pt($self,pt,m1,m2);
  }
  


  
/*####
in centroid.h
######*/
  
  /*calculates the centroid structure for alignment, and the distance to it*/
  double *centroid(char *OUTPUT)
  {
	  double *dist= malloc(sizeof(double));
	  *OUTPUT = vrna_centroid($self,dist);
	  return dist;
  }



  /*##########
   from constraints.h
################*/
  
  void contraints_add(const char *constraint, unsigned int options)
  {
	  vrna_constraints_add($self,constraint, options);
  }
  
  void hc_init()
  {
	  vrna_hc_init($self);
  }
  
  void hc_add_bp_nonspecific(int i,int d,char option)
  {
	vrna_hc_add_bp_nonspecific($self,i,d,option);  
  }
  
  void sc_init(vrna_fold_compound_t *vc)
  {
	  vrna_sc_init($self);
  }
  
  
  void sc_remove(){
    vrna_sc_remove($self);
  }

  void sc_add_up(const double *constraints, unsigned int options=VRNA_OPTION_MFE){
    vrna_sc_add_up($self, constraints, options);
  }

  void sc_add_bp(const double **constraints, unsigned int options=VRNA_OPTION_MFE){
    vrna_sc_add_bp($self, constraints, options);
  }

  int sc_add_hi_motif(const char *seq,const char *structure,double energy,unsigned int options=VRNA_OPTION_MFE){

    return vrna_sc_add_hi_motif($self, seq, structure, energy, options);
  }
  
  void sc_remove()
  {
	  vrna_sc_remove($self);
  }
  
  
  int sc_add_SHAPE_deigan(const double *reactivities,
                              double m,
                              double b,
                              unsigned int options)
  {
	  return vrna_sc_add_SHAPE_deigan($self,reactivities,m,b,options);
  }


  int sc_add_SHAPE_deigan_ali(const char **shape_files,
                                  const int *shape_file_association,
                                  double m,
                                  double b,
                                  unsigned int options)
  {
	  return vrna_sc_add_SHAPE_deigan_ali($self,shape_files,shape_file_association,m,b,options);                                  
  }

  int sc_add_SHAPE_zarringhalam(const double *reactivities,
                                    double b,
                                    double default_value,
                                    const char *shape_conversion,
                                    unsigned int options)
  {
	  return vrna_sc_add_SHAPE_zarringhalam($self,reactivities,b,default_value,shape_conversion,options);
  }
  
/*END constraints.h
#########################################*/
	
/*#######start ligand.h###*/

/*only double is argument energy of function, not float or double*/

int sc_add_hi_motif(const char *seq,
                      const char *structure,
                      double energy,
                      unsigned int options)
{
	return vrna_sc_add_hi_motif($self,seq,structure,energy,options);
}

/* return 1 or 0 if success and the base positions in a given structure if the given motif was found*/
%apply int *OUTPUT {int *i, int *j, int *k, int *l};  /* HERE more return parameters are defined*/
int sc_detect_hi_motif(const char *structure,
                        int *i,
                        int *j,
                        int *k,
                        int *l)
{
	return vrna_sc_detect_hi_motif($self,structure,i,j,k,l);
}

/* return 1 or 0 if success and the base positions in a givens tructure if the given motif was found*/
%apply int *OUTPUT {int *i, int *j, int *k, int *l};  /* HERE more return parameters are defined*/
int sc_get_hi_motif(    int *i,
                        int *j,
                        int *k,
                        int *l)
{
	
	return vrna_sc_get_hi_motif($self,i,j,k,l);
}







}


/*
 *  Rename all the preprocessor macros defined in data_structures.h
 *  (wrapped as constants)
 */
%constant unsigned char STATUS_MFE_PRE  = VRNA_STATUS_MFE_PRE;
%constant unsigned char STATUS_MFE_POST = VRNA_STATUS_MFE_POST;
%constant unsigned char STATUS_PF_PRE   = VRNA_STATUS_PF_PRE;
%constant unsigned char STATUS_PF_POST  = VRNA_STATUS_PF_POST;

%constant unsigned int OPTION_DEFAULT   = VRNA_OPTION_DEFAULT;
%constant unsigned int OPTION_MFE       = VRNA_OPTION_MFE;
%constant unsigned int OPTION_PF        = VRNA_OPTION_PF;
%constant unsigned int OPTION_HYBRID    = VRNA_OPTION_HYBRID;
%constant unsigned int OPTION_EVAL_ONLY = VRNA_OPTION_EVAL_ONLY;
%constant unsigned int OPTION_WINDOW    = VRNA_OPTION_WINDOW;

%rename(basepair) vrna_basepair_t;

typedef struct {
  int i;
  int j;
} vrna_basepair_t;

%include <ViennaRNA/data_structures.h>
