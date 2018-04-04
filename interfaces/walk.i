/****************************************************/
/* BEGIN interface for energy landscape exploration */
/****************************************************/

%extend vrna_fold_compound_t{

#include <vector>

#ifdef SWIGPYTHON
%feature("autodoc") path;
%feature("kwargs") path;
%feature("autodoc") path_gradient;
%feature("kwargs") path_gradient;
%feature("autodoc") path_random;
%feature("kwargs") path_random;
#endif

  std::vector<vrna_move_t>
  path(std::vector<int> &pt,
       unsigned int steps,
       unsigned int options = VRNA_PATH_DEFAULT)
  {
    int i;
    std::vector<short>::iterator it;
    std::vector<vrna_move_t>  v; /* fill vector with returned vrna_move_t */
    vrna_move_t *move_t, *ptr;
    std::vector<short> vc;

    transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);

    move_t = ptr = vrna_path($self, (short*)&vc[0], steps, options);
    
    if (ptr)
      while ((ptr->pos_5 != 0) && (ptr->pos_3 != 0)) {
        vrna_move_t m;
        m = vrna_move_init(ptr->pos_5, ptr->pos_3);
        v.push_back(m);
        ptr++;
      }

    /* copy over the values from vc to pt */
    for (i = 0, it = vc.begin(); it != vc.end(); ++it, i++)
      pt[i] = *it;

    free(move_t);
    return v;
  }

  std::vector<vrna_move_t>
  path_gradient(std::vector<int> &pt,
                unsigned int options = VRNA_PATH_DEFAULT)
  {
    int i;
    std::vector<short>::iterator it;
    std::vector<vrna_move_t>  v; /* fill vector with returned vrna_move_t */
    vrna_move_t *move_t, *ptr;

    std::vector<short> vc;
    transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);

    move_t = ptr = vrna_path_gradient($self, (short*)&vc[0], options);

    if (ptr)
      while ((ptr->pos_5 != 0) && (ptr->pos_3 != 0)) {
        vrna_move_t m;
        m = vrna_move_init(ptr->pos_5, ptr->pos_3);
        v.push_back(m);
        ptr++;
      }

    /* copy over the values from vc to pt */
    for (i = 0, it = vc.begin(); it != vc.end(); ++it, i++)
      pt[i] = *it;

    free(move_t);
    return v;
  }

  std::vector<vrna_move_t>
  path_random(std::vector<int> &pt,
              unsigned int steps,
              unsigned int options = VRNA_PATH_DEFAULT)
  {
    int i;
    std::vector<short>::iterator it;
    std::vector<vrna_move_t>  v; /* fill vector with returned vrna_move_t */
    vrna_move_t *move_t, *ptr;

    std::vector<short> vc;
    transform(pt.begin(), pt.end(), back_inserter(vc), convert_vecint2vecshort);

    move_t = ptr = vrna_path_random($self, (short*)&vc[0], steps, options);

    if (ptr)
      while ((ptr->pos_5 != 0) && (ptr->pos_3 != 0)) {
        vrna_move_t m;
        m = vrna_move_init(ptr->pos_5, ptr->pos_3);
        v.push_back(m);
        ptr++;
      }

    /* copy over the values from vc to pt */
    for (i = 0, it = vc.begin(); it != vc.end(); ++it, i++)
      pt[i] = *it;

    free(move_t);
    return v;
  }

}

%constant unsigned int PATH_STEEPEST_DESCENT      = VRNA_PATH_STEEPEST_DESCENT;
%constant unsigned int PATH_RANDOM                = VRNA_PATH_RANDOM;
%constant unsigned int PATH_NO_TRANSITION_OUTPUT  = VRNA_PATH_NO_TRANSITION_OUTPUT;
%constant unsigned int PATH_DEFAULT               = VRNA_PATH_DEFAULT;

%include <ViennaRNA/walk.h>
