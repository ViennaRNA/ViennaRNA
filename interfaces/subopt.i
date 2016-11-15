/**********************************************/
/* BEGIN interface for Suboptimal Structure   */
/* predictionn                                */
/**********************************************/

// from subopt.h


extern  int subopt_sorted;                       /* sort output by energy */

%nodefaultdtor SOLUTION;

typedef struct {
  float energy;
  char  *structure;
} SOLUTION;

%extend SOLUTION {
        SOLUTION *get(int i) {
           return self+i;
        }

        int size() {
           SOLUTION *s;
           for (s=self; s->structure; s++);
           return (int)(s-self);
        }

        ~SOLUTION() {
           SOLUTION *s;
           for (s=$self; s->structure; s++) free(s->structure);
           free($self);
        }
}

%rename (subopt) my_subopt;

/* We now have two different modes of return values for the subopt() method.
 * The 'old' one returns a single object of type 'SOLUTION' whereas the
 * 'new' return value actually is a std::vector<SOLUTION>. The latter enables
 * us to let swig generate native lists/arrays in the target language.
 *
 * In order to not confuse swig with deletion of a SOLUTION* type and
 * a single SOLUTION object, we just define a new struct 'subopt_solution' with
 * the same properties as SOLUTION, add a template for a std::vector<T> and
 * tell swig that we have our own destructor for objects of this type
 */
%{

extern "C" {
  typedef struct {
    float energy;
    char *structure;
  } subopt_solution;
}

%}

%nodefaultdtor subopt_solution;

typedef struct {
    float energy;
    char *structure;
} subopt_solution;

%extend subopt_solution {

  ~subopt_solution() {
    free($self->structure);
    free($self);
  }
}

namespace std {
  %template(SuboptVector) std::vector<subopt_solution>;
};


%{

  SOLUTION *my_subopt(char *seq, char *constraint, int delta, FILE *fp){
    return subopt(seq, constraint, delta, fp);
  }

  std::vector<subopt_solution> my_subopt(char *seq, int delta, FILE *fp){
    std::vector<subopt_solution> ret;
    SOLUTION *sol = subopt(seq, NULL, delta, fp);
    for(int i = 0; sol[i].structure != NULL; i++){
      subopt_solution a;
      a.energy = sol[i].energy;
      a.structure = sol[i].structure;
      ret.push_back(a);
    }
    free(sol);
    /* The memory occupied by the individual structures will be free'd automatically
       by swig, when the vector is destroyed
    */
    return ret;
  }

  std::vector<subopt_solution> my_subopt(char *seq, int delta){
    std::vector<subopt_solution> ret;
    SOLUTION *sol = subopt(seq, NULL, delta, NULL);
    for(int i = 0; sol[i].structure != NULL; i++){
      subopt_solution a;
      a.energy = sol[i].energy;
      a.structure = sol[i].structure;
      ret.push_back(a);
    }
    free(sol);
    /* The memory occupied by the individual structures will be free'd automatically
       by swig, when the vector is destroyed
    */
    return ret;
  }

  SOLUTION *my_subopt(char *seq, char *constraint, int delta){
    return subopt(seq, constraint, delta, NULL);
  }

%}

%newobject my_subopt;

SOLUTION *my_subopt(char *seq, char *constraint, int delta, FILE *fp);
std::vector<subopt_solution> my_subopt(char *seq, int delta, FILE *fp);
std::vector<subopt_solution> my_subopt(char *seq, int delta);
SOLUTION *my_subopt(char *seq, char *constraint, int delta);

%extend vrna_fold_compound_t {

  std::vector<subopt_solution> subopt(int delta, int sorted=1){
    std::vector<subopt_solution> ret;
    SOLUTION *sol = vrna_subopt($self, delta, sorted, NULL);
    for(int i = 0; sol[i].structure != NULL; i++){
      subopt_solution a;
      a.energy = sol[i].energy;
      a.structure = sol[i].structure;
      ret.push_back(a);
    }
    free(sol);
    /* The memory occupied by the individual structures will be free'd automatically
       by swig, when the vector is destroyed
    */
    return ret;
  }

  std::vector<subopt_solution> subopt(int delta, int sorted, FILE *fp){
    std::vector<subopt_solution> ret;
    SOLUTION *sol = vrna_subopt($self, delta, sorted, fp);
    for(int i = 0; sol[i].structure != NULL; i++){
      subopt_solution a;
      a.energy = sol[i].energy;
      a.structure = sol[i].structure;
      ret.push_back(a);
    }
    free(sol);
    /* The memory occupied by the individual structures will be free'd automatically
       by swig, when the vector is destroyed
    */
    return ret;
  }

  std::vector<subopt_solution> subopt_zuker(void){
    std::vector<subopt_solution> ret;
    SOLUTION *sol = vrna_subopt_zuker($self);
    for(int i = 0; sol[i].structure != NULL; i++){
      subopt_solution a;
      a.energy = sol[i].energy;
      a.structure = sol[i].structure;
      ret.push_back(a);
    }
    free(sol);
    /* The memory occupied by the individual structures will be free'd automatically
       by swig, when the vector is destroyed
    */
    return ret;
  }
}

%ignore subopt;
%ignore subopt_par;
%ignore subopt_circ;
/*
%ignore zukersubopt;
*/
%ignore zukersubopt_par;

%include  <ViennaRNA/subopt.h>
