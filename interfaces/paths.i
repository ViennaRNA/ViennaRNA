/**********************************************/
/* BEGIN interface for (re-)folding path      */
/* implementations                            */
/**********************************************/

#include <iostream>
#include <vector>
#include <limits.h>

%ignore path;

/* scripting language access through 'fold_compound' instead of 'vrna_fold_compound_t' */
%rename(path) vrna_path_t;

/* no default constructor / destructor */
%nodefaultdtor vrna_path_t;

typedef struct {
  unsigned int type;
  double en;  /**<  @brief  Free energy of current structure */
  char *s;    /**<  @brief  Secondary structure in dot-bracket notation */
  vrna_move_t   move;
} vrna_path_t;

%extend vrna_path_t {
  ~vrna_path_t() {
    free($self->s);
    free($self);
  }
}

/**********************************************/

%ignore vrna_path_options_t;
%ignore vrna_path_options_s;

%rename(path_options) vrna_path_options_s;

typedef struct {} vrna_path_options_s;

%nodefaultdtor vrna_path_options_s;
%nodefaultctor vrna_path_options_s;

%extend vrna_path_options_s {
  vrna_path_options_s()
  {
    return NULL;
  }

  ~vrna_path_options_s(){
    vrna_path_options_free($self);
  }
}

/**********************************************/

/* Provide an overloaded wrapper for vrna_path_options_findpath() */
%rename(path_options_findpath) my_path_options_findpath;

%{
  struct vrna_path_options_s *
  my_path_options_findpath(int          width = 10,
                           unsigned int type  = VRNA_PATH_TYPE_DOT_BRACKET)
  {
    return vrna_path_options_findpath(width, type);
  }
%}

struct vrna_path_options_s *
my_path_options_findpath(int          width = 10,
                         unsigned int type  = VRNA_PATH_TYPE_DOT_BRACKET);

#ifdef SWIGPYTHON
%feature("autodoc") my_path_options_findpath;
%feature("kwargs") my_path_options_findpath;
#endif

%ignore vrna_path_options_findpath;

/**********************************************/

/* Provide a wrapper for get_path() */
%rename (get_path) my_get_path;

%{
  std::vector<vrna_path_t> my_get_path(std::string seq,
                                       std::string s1,
                                       std::string s2,
                                       int         maxkeep)
  {
    std::vector<vrna_path_t>  v; /* fill vector with returned vrna_path_t*/
    vrna_path_t *path_t, *ptr;

    path_t = ptr = get_path(seq.c_str(), s1.c_str(), s2.c_str(), maxkeep);

    while (ptr->s != NULL)
    {
        vrna_path_t p;
        p.en = ptr->en;
        p.s  = ptr->s;
        v.push_back(p);
        ptr++;
        
    }
    free(path_t);
    return v;
  }
%}

#ifdef SWIGPYTHON
%feature("autodoc") my_get_path;
%feature("kwargs") my_get_path;
#endif

std::vector<vrna_path_t> my_get_path(std::string seq, std::string s1, std::string s2, int maxkeep);
%ignore get_path;

/**********************************************/

%extend vrna_fold_compound_t{
  
#ifdef SWIGPYTHON
%feature("autodoc") path_findpath;
%feature("kwargs") path_findpath;
%feature("autodoc") path_findpath_saddle;
%feature("kwargs") path_findpath_saddle;

  PyObject *
  path_findpath_saddle(std::string s1, std::string s2, int width = 1, int maxE = INT_MAX){
    PyObject *E_obj = Py_None;

    int E = vrna_path_findpath_saddle_ub($self, s1.c_str(), s2.c_str(), width, maxE);

    if (E < maxE)
      E_obj = Py_BuildValue("i", E);
    else
      Py_INCREF(Py_None); /* increase reference count for Py_None */

    return E_obj;
  }
#endif

#ifdef SWIGPERL5
  SV *
  path_findpath_saddle(std::string s1, std::string s2, int width = 1, int maxE = INT_MAX - 1){
    SV *E_obj;

    int E = vrna_path_findpath_saddle_ub($self, s1.c_str(), s2.c_str(), width, maxE);

    if (E < maxE)
      E_obj = newSViv((IV)E);
    else
      E_obj = newSV(0);

    sv_2mortal(E_obj);

    return E_obj;
  }

#endif

  std::vector<vrna_path_t> path_findpath(std::string s1, std::string s2, int width = 1, int maxE = INT_MAX - 1){
      std::vector<vrna_path_t>  v; /* fill vector with returned vrna_path_t*/
      vrna_path_t *path_t, *ptr;
      path_t = ptr = vrna_path_findpath_ub($self, s1.c_str(), s2.c_str(), width, maxE);

      if (ptr) {
        while (ptr->s != NULL)
        {
            vrna_path_t p;
            p.en = ptr->en;
            p.s  = ptr->s;
            v.push_back(p);
            ptr++;
        }
        free(path_t);
      }
      return v;
  }

#ifdef SWIGPYTHON
%feature("autodoc") path_direct;
%feature("kwargs") path_direct;
#endif

  std::vector<vrna_path_t>
  path_direct(std::string                 s1,
              std::string                 s2,
              int                         maxE = INT_MAX - 1,
              struct vrna_path_options_s  *options = NULL)
  {
      std::vector<vrna_path_t>  v; /* fill vector with returned vrna_path_t*/
      vrna_path_t *path_t, *ptr;
      path_t = ptr = vrna_path_direct_ub($self, s1.c_str(), s2.c_str(), maxE, options);

      if (ptr) {
        if (ptr->type == VRNA_PATH_TYPE_DOT_BRACKET)
          for (; ptr->s != NULL; ptr++) {
            vrna_path_t p;
            p.type = ptr->type;
            p.en   = ptr->en;
            p.s    = ptr->s;
            p.move = ptr->move;
            v.push_back(p);
          }
        else if (ptr->type == VRNA_PATH_TYPE_MOVES)
          for (; ptr->move.pos_5 != 0; ptr++) {
            vrna_path_t p;
            p.type = ptr->type;
            p.en   = ptr->en;
            p.s    = ptr->s;
            p.move = ptr->move;
            v.push_back(p);
          }
      }

      free(path_t);

      return v;
  }

}

/**********************************************/

%constant unsigned int PATH_TYPE_DOT_BRACKET  = VRNA_PATH_TYPE_DOT_BRACKET;
%constant unsigned int PATH_TYPE_MOVES        = VRNA_PATH_TYPE_MOVES;

%include <ViennaRNA/landscape/paths.h>
%include <ViennaRNA/landscape/findpath.h>

