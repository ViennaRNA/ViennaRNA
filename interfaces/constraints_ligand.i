/**********************************************/
/* BEGIN interface for structure ligand constraints */
/**********************************************/

%extend vrna_fold_compound_t {

int sc_add_hi_motif(const char *seq,
                    const char *structure,
                    FLT_OR_DBL energy,
                    unsigned int options=VRNA_OPTION_DEFAULT){

  return vrna_sc_add_hi_motif($self,seq,structure,energy,options);
}


/*
%apply int *OUTPUT {int *i, int *j, int *k, int *l};  /* HERE more return parameters are defined
int sc_detect_hi_motif(const char *structure,
                        int *i,
                        int *j,
                        int *k,
                        int *l)*/
//   int sc_detect_hi_motif(const char *structure)                      
// {
//         
//         
//         int ret = vrna_sc_detect_hi_motif($self,structure,&i,&j,&k,&l);
//         
//         return ret;
// }



/* 
// %apply int *OUTPUT {int *i, int *j, int *k, int *l};  /* HERE more return parameters are defined*/
// int sc_get_hi_motif(    int *i,
//                         int *j,
//                         int *k,
//                         int *l)
// {
//         
//         int ret =  vrna_sc_get_hi_motif($self,i,j,k,l);
//         return ret;
// }

}

%include  <ViennaRNA/constraints_ligand.h>
