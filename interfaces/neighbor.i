/****************************************************/
/* BEGIN interface for neighbor generation          */
/****************************************************/

#ifdef SWIGPERL5
%rename(_next) vrna_move_s::next;
#endif

/* scripting language access through 'move' instead of 'vrna_move_t' */
%rename(move) vrna_move_t;

/* no default destructor */
%nodefaultctor vrna_move_t;
%nodefaultdtor vrna_move_t;

typedef struct {
  int pos_5;
  int pos_3;
} vrna_move_t;


/* create object oriented interface for vrna_move_t */
%extend vrna_move_t {
  vrna_move_t() {
    vrna_move_t *m = (vrna_move_t *)vrna_alloc(sizeof(vrna_move_t));
    *m = vrna_move_init(0, 0);
    return m;
  }

  vrna_move_t(int pos_5, int pos_3) {
    vrna_move_t *m = (vrna_move_t *)vrna_alloc(sizeof(vrna_move_t));
    *m = vrna_move_init(pos_5, pos_3);
    return m;
  }

  ~vrna_move_t(){
    vrna_move_list_free($self->next);
    free($self);
  }
}

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


%constant unsigned int MOVESET_INSERTION  = VRNA_MOVESET_INSERTION;
%constant unsigned int MOVESET_DELETION   = VRNA_MOVESET_DELETION;
%constant unsigned int MOVESET_SHIFT      = VRNA_MOVESET_SHIFT;
%constant unsigned int MOVESET_NO_LP      = VRNA_MOVESET_NO_LP;
%constant unsigned int MOVESET_DEFAULT    = VRNA_MOVESET_DEFAULT;

%include <ViennaRNA/neighbor.h>
