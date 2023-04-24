/****************************************************/
/* BEGIN interface for neighbor generation          */
/****************************************************/

%apply  std::vector<vrna_move_t> *OUTPUT { std::vector<vrna_move_t> *invalid_moves};

%extend vrna_fold_compound_t{

#ifdef SWIGPYTHON
%feature("autodoc") neighbors;
%feature("kwargs") neighbors;

  std::vector<vrna_move_t>
  neighbors(var_array<short> &pt,
            unsigned int      options = VRNA_MOVESET_DEFAULT)
  {
    std::vector<vrna_move_t>  v; /* fill vector with returned vrna_move_t */
    vrna_move_t *move_t, *ptr;

    move_t = ptr = vrna_neighbors($self, pt.data, options);

    if (ptr)
      while ((ptr->pos_5 != 0) && (ptr->pos_3 != 0)) {
        vrna_move_t m;
        m = vrna_move_init(ptr->pos_5, ptr->pos_3);
        v.push_back(m);
        ptr++;
      }

    free(move_t);
    return v;
  }

  %newobject vrna_fold_compound_t::move_neighbor_diff;

  std::vector<vrna_move_t>
  move_neighbor_diff(var_array<short> &pt,
                     vrna_move_t      move,
                     std::vector<vrna_move_t> *invalid_moves,
                     unsigned int          options = VRNA_MOVESET_DEFAULT)
  {
    std::vector<vrna_move_t>  v; /* fill vector with returned vrna_move_t */
    vrna_move_t *move_t, *ptr, *ptr2, *mv_invalid;

    mv_invalid = NULL;
    move_t = ptr = vrna_move_neighbor_diff($self,
                                           pt.data,
                                           move,
                                           &mv_invalid,
                                           options);

    if (mv_invalid) {
      *invalid_moves = std::vector<vrna_move_t>();
      ptr2 = mv_invalid;
      while ((ptr2->pos_5 != 0) && (ptr2->pos_3 != 0)) {
        vrna_move_t m;
        m = vrna_move_init(ptr2->pos_5, ptr2->pos_3);
        (*invalid_moves).push_back(m);
        ptr2++;
      }
    }

    if (ptr) {
      while ((ptr->pos_5 != 0) && (ptr->pos_3 != 0)) {
        vrna_move_t m;
        m = vrna_move_init(ptr->pos_5, ptr->pos_3);
        v.push_back(m);
        ptr++;
      }
    }

    free(move_t);
    free(mv_invalid);

    return v;
  }

}

#else
  std::vector<vrna_move_t>
  neighbors(std::vector<int>  pt,
            unsigned int      options = VRNA_MOVESET_DEFAULT)
  {
    std::vector<vrna_move_t>  v; /* fill vector with returned vrna_move_t */
    vrna_move_t *move_t, *ptr;
    std::vector<short> vc;

    transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);

    move_t = ptr = vrna_neighbors($self, (short*)&vc[0], options);

    if (ptr)
      while ((ptr->pos_5 != 0) && (ptr->pos_3 != 0)) {
        vrna_move_t m;
        m = vrna_move_init(ptr->pos_5, ptr->pos_3);
        v.push_back(m);
        ptr++;
      }

    free(move_t);
    return v;
  }
}

#endif

%constant unsigned int NEIGHBOR_CHANGE  = VRNA_NEIGHBOR_CHANGE;
%constant unsigned int NEIGHBOR_INVALID = VRNA_NEIGHBOR_INVALID;
%constant unsigned int NEIGHBOR_NEW     = VRNA_NEIGHBOR_NEW;


%include <ViennaRNA/landscape/neighbor.h>

