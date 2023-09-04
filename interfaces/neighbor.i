/****************************************************/
/* BEGIN interface for neighbor generation          */
/****************************************************/


%newobject vrna_fold_compound_t::neighbors;
%newobject vrna_fold_compound_t::move_neighbor_diff;

%extend vrna_fold_compound_t{

#ifdef SWIGPYTHON
%apply  var_array< vrna_move_t > *&OUTPUT { var_array< vrna_move_t > *&invalid_moves };

%feature("autodoc") neighbors;
%feature("kwargs") neighbors;


  var_array<vrna_move_t> *
  neighbors(var_array<short> &pt,
            unsigned int     options = VRNA_MOVESET_DEFAULT)
  {
    var_array<vrna_move_t>  *v;
    vrna_move_t             *move_t, *ptr;

    v      = NULL;
    move_t = vrna_neighbors($self, pt.data, options);

    if (move_t) {
      /* get size of new and changed moves */
      size_t n = 0;
      for (ptr = move_t; ptr->pos_5 != 0; ptr++, n++);
      v = var_array_new(n, move_t, VAR_ARRAY_LINEAR | VAR_ARRAY_OWNED);
    }

    return v;
  }


  var_array<vrna_move_t> *
  move_neighbor_diff(var_array<short> &pt,
                     vrna_move_t      move,
                     var_array<vrna_move_t> *&invalid_moves,
                     unsigned int          options = VRNA_MOVESET_DEFAULT)
  {
    var_array<vrna_move_t> *v; /* fill vector with returned vrna_move_t */
    vrna_move_t *move_t, *mv_invalid, *ptr;

    mv_invalid  = NULL;
    v           = NULL;
    move_t = vrna_move_neighbor_diff($self,
                                     pt.data,
                                     move,
                                     &mv_invalid,
                                     options);

    if (mv_invalid) {
      /* get size of invalid moves */
      size_t n = 0;
      for (ptr = mv_invalid; ptr->pos_5 != 0; ptr++, n++);
      invalid_moves = var_array_new(n, mv_invalid, VAR_ARRAY_LINEAR | VAR_ARRAY_OWNED);
    }

    if (move_t) {
      /* get size of new and changed moves */
      size_t n = 0;
      for (ptr = move_t; ptr->pos_5 != 0; ptr++, n++);
      v = var_array_new(n, move_t, VAR_ARRAY_LINEAR | VAR_ARRAY_OWNED);
    }

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

