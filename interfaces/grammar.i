
/**********************************************/
/* BEGIN interface for grammar extensions     */
/**********************************************/

%extend vrna_fold_compound_t {
#include <string>
#include <cstring>

  void ud_add_motif(std::string motif, double motif_en, unsigned int options=VRNA_UNSTRUCTURED_DOMAIN_ALL_LOOPS){
    vrna_ud_add_motif($self, motif.c_str(), motif_en, options);
  }

  void ud_remove(void){
    vrna_ud_remove($self);
  }
}

%constant unsigned int UNSTRUCTURED_DOMAIN_EXT_LOOP   = VRNA_UNSTRUCTURED_DOMAIN_EXT_LOOP;
%constant unsigned int UNSTRUCTURED_DOMAIN_HP_LOOP    = VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP;
%constant unsigned int UNSTRUCTURED_DOMAIN_INT_LOOP   = VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP;
%constant unsigned int UNSTRUCTURED_DOMAIN_MB_LOOP    = VRNA_UNSTRUCTURED_DOMAIN_MB_LOOP;
%constant unsigned int UNSTRUCTURED_DOMAIN_ALL_LOOPS  = VRNA_UNSTRUCTURED_DOMAIN_ALL_LOOPS;
%constant unsigned int UNSTRUCTURED_DOMAIN_MOTIF      = VRNA_UNSTRUCTURED_DOMAIN_MOTIF;

%include  <ViennaRNA/grammar.h>
%include  <ViennaRNA/unstructured_domains.h>
%include  <ViennaRNA/structured_domains.h>
