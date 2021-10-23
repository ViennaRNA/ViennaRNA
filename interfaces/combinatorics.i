/******************************************************/
/* BEGIN interface for Combinatorics Implementations  */
/******************************************************/

%ignore vrna_enumerate_necklaces;
%ignore vrna_rotational_symmetry_num;
%ignore vrna_rotational_symmetry;

%rename (enumerate_necklaces) my_enumerate_necklaces;
%rename (rotational_symmetry) my_rotational_symmetry;

%{
  std::vector<std::vector<int> >
  my_enumerate_necklaces( std::vector<unsigned int> entity_counts)
  {
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


  std::vector<unsigned int>
  my_rotational_symmetry(std::string string)
  {
    std::vector<unsigned int> positions;
    unsigned int i, r, *pos;
    
    r = vrna_rotational_symmetry_pos(string.c_str(), &pos);

    if (r)
      for (i = 0; i < r; i++)
        positions.push_back(pos[i]);

    free(pos);

    return positions;
  }


  std::vector<unsigned int>
  my_rotational_symmetry(std::vector<unsigned int> string)
  {
    std::vector<unsigned int> positions;
    unsigned int i, r, *pos;
    
    r = vrna_rotational_symmetry_pos_num((unsigned int*)&string[0], string.size(), &pos);

    if (r)
      for (i = 0; i < r; i++)
        positions.push_back(pos[i]);

    free(pos);

    return positions;
  }

%}

std::vector<std::vector<int> > my_enumerate_necklaces( std::vector<unsigned int> entity_counts);
std::vector<unsigned int>      my_rotational_symmetry(std::vector<unsigned int> string);
std::vector<unsigned int>      my_rotational_symmetry(std::string string);


%extend vrna_fold_compound_t {

  std::vector<unsigned int>
  rotational_symmetry_db(std::string structure) {
    std::vector<unsigned int> positions;
    unsigned int i, r, *pos;

    r = vrna_rotational_symmetry_db_pos($self, structure.c_str(), &pos);

    if (r)
      for (i = 0; i < r; i++)
        positions.push_back(pos[i]);

    free(pos);

    return positions;
  }
}


%{
  std::vector<unsigned int>
  boustrophedon(unsigned int start,
                unsigned int end)
  {
    std::vector<unsigned int> sequence;
    unsigned int *result = vrna_boustrophedon(start, end);

    if (result) {
      for (size_t i = 0; i <= result[0]; i++)
        sequence.push_back(result[i]);

      free(result);
    }

    return sequence;
  }

  unsigned int
  boustrophedon(unsigned int start,
                unsigned int end,
                unsigned int pos)
  {
    return vrna_boustrophedon_pos((size_t)start,
                                  (size_t)end,
                                  (size_t)pos);
  }
%}

std::vector<unsigned int> boustrophedon(unsigned int start, unsigned int end);
unsigned int              boustrophedon(unsigned int start, unsigned int end, unsigned int pos);


%include  <ViennaRNA/combinatorics.h>
