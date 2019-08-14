/****************************************************/
/* BEGIN interface for neighbor generation          */
/****************************************************/

%extend vrna_fold_compound_t{

#ifdef SWIGPYTHON
%feature("autodoc") neighbors;
%feature("kwargs") neighbors;
#endif

  std::vector<vrna_move_t>
  neighbors(std::vector<int> pt,
            unsigned int options = VRNA_MOVESET_DEFAULT)
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

%constant unsigned int NEIGHBOR_CHANGE  = VRNA_NEIGHBOR_CHANGE;
%constant unsigned int NEIGHBOR_INVALID = VRNA_NEIGHBOR_INVALID;
%constant unsigned int NEIGHBOR_NEW     = VRNA_NEIGHBOR_NEW;


%include <ViennaRNA/landscape/neighbor.h>

