/******************************************************/
/* BEGIN interface for Combinatorics Implementations  */
/******************************************************/

%ignore vrna_enumerate_necklaces;

%rename (enumerate_necklaces) my_enumerate_necklaces;

%{
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
%}

std::vector<std::vector<int> > my_enumerate_necklaces( std::vector<unsigned int> entity_counts);

%include  <ViennaRNA/combinatorics.h>
