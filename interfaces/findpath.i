/**********************************************/
/* BEGIN interface for findpath heursitic     */
/**********************************************/

#include <iostream>
#include <vector>

/* scripting language access through 'fold_compound' instead of 'vrna_fold_compound_t' */
%rename(path) vrna_path_t;

/* no default constructor / destructor */
%nodefaultdtor vrna_path_t;

typedef struct {
  double en;  /**<  @brief  Free energy of current structure */
  char *s;    /**<  @brief  Secondary structure in dot-bracket notation */
} vrna_path_t;

%extend vrna_path_t {
  vrna_path_t *get(int i) {
    return $self+i;
  }

  int size() {
    path_t *st;
    for (st=$self; st->s; st++);
    return (int)(st-$self);
  }

  ~vrna_path_t() {
    vrna_path_t *st;
    for (st=$self; st->s; st++)
      free(st->s);
    free($self);
  }
}


%extend vrna_fold_compound_t{
  
  int path_findpath_saddle(const char *struc1, const char *struc2, int max){
    return vrna_path_findpath_saddle($self,struc1,struc2,max);
  }
  
  std::vector<vrna_path_t> path_findpath(std::string s1, std::string s2, int maxkeep){
    
      std::vector<vrna_path_t>  v; /* fill vector with returned vrna_path_t*/
      vrna_path_t *path_t, *ptr;
      path_t = ptr = vrna_path_findpath($self,s1.c_str(),s2.c_str(),maxkeep);
      
      while(ptr->s != NULL)
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

