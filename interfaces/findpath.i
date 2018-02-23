/**********************************************/
/* BEGIN interface for findpath heursitic     */
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
  double en;  /**<  @brief  Free energy of current structure */
  char *s;    /**<  @brief  Secondary structure in dot-bracket notation */
} vrna_path_t;

%extend vrna_path_t {
  ~vrna_path_t() {
    free($self->s);
    free($self);
  }
}

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
std::vector<vrna_path_t> my_get_path(std::string seq, std::string s1, std::string s2, int maxkeep);
%ignore get_path;


%extend vrna_fold_compound_t{
  
#ifdef SWIGPYTHON
%feature("kwargs") path_findpath_saddle;

  PyObject *
  path_findpath_saddle(std::string s1, std::string s2, int width = 1, int maxE = INT_MAX){
    PyObject *E_obj = Py_None;

    int E = vrna_path_findpath_saddle_ub($self, s1.c_str(), s2.c_str(), width, maxE);

    if (E < maxE)
      E_obj = Py_BuildValue("i", E);

    return E_obj;
  }
#endif

#ifdef SWIGPERL5
  SV *
  path_findpath_saddle(std::string s1, std::string s2, int width = 1, int maxE = INT_MAX){
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

  std::vector<vrna_path_t> path_findpath(std::string s1, std::string s2, int width = 1, int maxE = INT_MAX){
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

}


%include <ViennaRNA/findpath.h>;

