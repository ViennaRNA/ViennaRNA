/*
    sequence.c
    
    Code for handling nucleotide sequences
    
    Part of the ViennaRNA Package
    (c) 2016 Ronny Lorenz
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/utils.h"
#include "ViennaRNA/string_utils.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/sequence.h"

/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE void set_sequence(vrna_seq_t *obj,
                          const char *string,
                          unsigned int options);

PRIVATE void free_sequence_data( vrna_seq_t *obj);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/
PUBLIC vrna_seq_t *
vrna_sequence(const char *string,
              unsigned int options){

  vrna_seq_t *data = NULL;

  if(string){
    data = (vrna_seq_t *) vrna_alloc(sizeof(vrna_seq_t));
    set_sequence(data, string, options);
  }

  return data;
}


PUBLIC int
vrna_sequence_add(vrna_fold_compound_t *vc,
                  const char *string,
                  unsigned int options){

  int ret = 0;

  if(vc && string){
    vc->nucleotides = (vrna_seq_t *) vrna_realloc(vc->nucleotides, sizeof(vrna_seq_t) * (vc->strands + 1));
    set_sequence(&(vc->nucleotides[vc->strands]), string, options);
    vc->strands++;

    ret = 1;
  }

  return ret;
}


PUBLIC int
vrna_sequence_remove( vrna_fold_compound_t *vc,
                      unsigned int i){

  unsigned int j, size;
  int ret = 0;

  if(vc){
    if(i < vc->strands){
      free_sequence_data(&(vc->nucleotides[i]));
      /* roll all nucleotide sequences behind the deleted one to the front */
      size = vc->strands - i - 1;
      if(size > 0)
        memmove(vc->nucleotides + i, vc->nucleotides + i + 1, sizeof(vrna_seq_t) * size);

      vc->strands--;
      vc->nucleotides = (vrna_seq_t *) vrna_realloc(vc->nucleotides, sizeof(vrna_seq_t) * vc->strands);

      ret = 1;
    }
  }

  return ret;
}


PUBLIC void
vrna_sequence_remove_all( vrna_fold_compound_t *vc){

  unsigned int i;

  if(vc){
    for(i = 0; i < vc->strands; i++)
      free_sequence_data(&(vc->nucleotides[i]));

    free(vc->nucleotides);

    vc->nucleotides = NULL;
    vc->strands     = 0;
  }
}

PRIVATE void
set_sequence( vrna_seq_t *obj,
              const char *string,
              unsigned int options){

  obj->string = strdup(string);
  vrna_seq_toupper(obj->string);
  obj->length = strlen(obj->string);

  switch(options){
    default:  obj->type = VRNA_SEQ_RNA;
  }

  obj->encoding = vrna_seq_encode(obj->string, NULL);
}

PRIVATE void
free_sequence_data( vrna_seq_t *obj){

  free(obj->string);
  free(obj->encoding);
  obj->type   = VRNA_SEQ_UNKNOWN;
  obj->length = 0;
}
