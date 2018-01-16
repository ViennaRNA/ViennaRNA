/**********************************************/
/* BEGIN interface for findpath heursitic     */
/**********************************************/

#include <iostream>
#include <vector>

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
  
  int path_findpath_saddle(const char *struc1, const char *struc2, int max){
    return vrna_path_findpath_saddle($self,struc1,struc2,max);
  }
  
  std::vector<vrna_path_t> path_findpath(std::string s1, std::string s2, int maxkeep){
    
      std::vector<vrna_path_t>  v; /* fill vector with returned vrna_path_t*/
      vrna_path_t *path_t, *ptr;
      path_t = ptr = vrna_path_findpath($self,s1.c_str(),s2.c_str(),maxkeep);
      
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

}


%include <ViennaRNA/findpath.h>;

