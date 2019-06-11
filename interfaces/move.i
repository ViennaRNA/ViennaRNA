/****************************************************/
/* BEGIN interface for structure moves              */
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

  int
  is_removal()
  {
    return vrna_move_is_removal((const vrna_move_t *)$self);
  }

  int
  is_insertion()
  {
    return vrna_move_is_insertion((const vrna_move_t *)$self);
  }
  
  int
  is_shift()
  {
    return vrna_move_is_shift((const vrna_move_t *)$self);
  }
  
  int
  compare(const vrna_move_t       *b,
          const std::vector<int>  pt = std::vector<int>())
  {
    int result;
    std::vector<short> vs;
    transform(pt.begin(), pt.end(), back_inserter(vs), convert_vecint2vecshort);

    result =  vrna_move_compare($self,
                                b,
                                (const short *)&vs[0]);

    return result;
  }
}

%constant unsigned int MOVESET_INSERTION  = VRNA_MOVESET_INSERTION;
%constant unsigned int MOVESET_DELETION   = VRNA_MOVESET_DELETION;
%constant unsigned int MOVESET_SHIFT      = VRNA_MOVESET_SHIFT;
%constant unsigned int MOVESET_NO_LP      = VRNA_MOVESET_NO_LP;
%constant unsigned int MOVESET_DEFAULT    = VRNA_MOVESET_DEFAULT;

%include  <ViennaRNA/landscape/move.h>
