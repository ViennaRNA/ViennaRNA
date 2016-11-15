/******************************************************/
/* BEGIN interface for Combinatorics Implementations  */
/******************************************************/

%ignore vrna_enumerate_necklaces;

%rename (enumerate_necklaces) my_enumerate_necklaces;

%{
#ifdef SWIGPERL
  SV *my_enumerate_necklaces( std::vector<unsigned int> entity_counts){
    /* add a 0 entry, just in case it has been forgotten */
    entity_counts.push_back(0);
    AV *myav  = NULL;
    SV *ret   = &PL_sv_undef;

    unsigned int **r = vrna_enumerate_necklaces((const unsigned int *)&entity_counts[0]);

    if(r){
      /* get line length */
      unsigned int n = 0;
      for(std::vector<unsigned int>::iterator it = entity_counts.begin(); it != entity_counts.end(); ++it)
        n += *it;

      SV **svs;
      int i, len;
      /* Figure out how many elements we have */
      for(len = 0; r[len]; len++);

      svs = (SV **) malloc(len*sizeof(SV *)); 
      for (i = 0; i < len ; i++) {
        /* create new array with current permutation */
        AV *perm = newAV();
        for(unsigned int j = 1; j <= n; j++)
          av_push(perm, newSVuv((UV) r[i][j]));

        /* store reference to current permutation in array */
        svs[i] = sv_2mortal(newRV_inc((SV*) perm));
        /* get rid of array ? */
      }
      myav =  av_make(len,svs);
      free(svs);
      free(r);
      ret = sv_2mortal(newRV_inc((SV*) myav));
    }

    return ret;
  }
#else
  std::vector<std::vector<int> >
  my_enumerate_necklaces( std::vector<unsigned int> entity_counts){
    std::vector<std::vector<int> > permutations;
    /* add a 0 entry, just in case it has been forgotten */
    entity_counts.push_back(0);
    unsigned int **result = vrna_enumerate_necklaces((const unsigned int *)&entity_counts[0]);
    if(result){
      /* get line length */
      unsigned int n = 0;
      for(std::vector<unsigned int>::iterator it = entity_counts.begin(); it != entity_counts.end(); ++it)
        n += *it;

      for(int i = 0; result[i]; i++){
        std::vector<int> line;
        for(unsigned int j = 1; j <= n; j++)
          line.push_back((int)result[i][j]);
        free(result[i]);
        permutations.push_back(line);
      }
      free(result);
    }
    return permutations;
  }
#endif
%}
#ifdef SWIGPERL
SV *my_enumerate_necklaces( std::vector<unsigned int> entity_counts);
#else
std::vector<std::vector<int> > my_enumerate_necklaces( std::vector<unsigned int> entity_counts);
#endif

%include  <ViennaRNA/combinatorics.h>
