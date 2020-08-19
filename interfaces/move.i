/****************************************************/
/* BEGIN interface for structure moves              */
/****************************************************/

%{
#include <sstream>
%}

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

#ifdef SWIGPYTHON
%feature("autodoc")vrna_move_t::vrna_move_t;
%feature("kwargs")vrna_move_t::vrna_move_t;
#endif

  vrna_move_t(int pos_5 = 0,
              int pos_3 = 0)
  {
    vrna_move_t *m = (vrna_move_t *)vrna_alloc(sizeof(vrna_move_t));
    *m = vrna_move_init(pos_5, pos_3);
    return m;
  }

  ~vrna_move_t()
  {
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

#ifdef SWIGPYTHON
  std::string
  __str__()
  {
    std::ostringstream out;
    out << "{ pos_5: " << $self->pos_5;
    out << ", pos_3: " << $self->pos_3;
    out << " }";

    return std::string(out.str());
  }

%pythoncode %{
def __repr__(self):
    # reformat string representation (self.__str__()) to something
    # that looks like a constructor argument list
    strthis = self.__str__().replace(": ", "=").replace("{ ", "").replace(" }", "")
    return  "%s.%s(%s)" % (self.__class__.__module__, self.__class__.__name__, strthis) 
%}
#endif

}

%constant unsigned int MOVESET_INSERTION  = VRNA_MOVESET_INSERTION;
%constant unsigned int MOVESET_DELETION   = VRNA_MOVESET_DELETION;
%constant unsigned int MOVESET_SHIFT      = VRNA_MOVESET_SHIFT;
%constant unsigned int MOVESET_NO_LP      = VRNA_MOVESET_NO_LP;
%constant unsigned int MOVESET_DEFAULT    = VRNA_MOVESET_DEFAULT;

%include  <ViennaRNA/landscape/move.h>
