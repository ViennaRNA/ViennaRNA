/**********************************************/
/* BEGIN interface for Suboptimal Structure   */
/* predictionn                                */
/**********************************************/

// from subopt.h


extern  int subopt_sorted;                       /* sort output by energy */

typedef struct {
  float energy;
  char  *structure;
} SOLUTION;

%extend SOLUTION {
        SOLUTION *get(int i) {
//           static int size=-1;
//           if (size<0) {
//             SOLUTION *s;
//             for (s=self; s->structure; s++);
//             size= (int) (s-self);
//           }
//           if (i>=size) {
//             warn("value out of range");
//             return NULL;
//           }
           return self+i;
        }

        int size() {
           SOLUTION *s;
           for (s=self; s->structure; s++);
           return (int)(s-self);
        }

        ~SOLUTION() {
           SOLUTION *s;
           for (s=self; s->structure; s++) free(s->structure);
           free(self);
        }
}

%rename (subopt) my_subopt;

%{
  SOLUTION *my_subopt(char *seq, char *constraint, int delta, FILE *fp){
    return subopt(seq, constraint, delta, fp);
  }
  std::vector<SOLUTION> my_subopt(char *seq, int delta, FILE *fp){
    std::vector<SOLUTION> ret;
    SOLUTION *sol = subopt(seq, NULL, delta, fp);
    for(int i = 0; sol[i].structure != NULL; i++){
      ret.push_back(sol[i]);
    }
    free(sol);
    /* hopefully the structures will be free'd when the vector is destroyed */
    return ret;
  }
  std::vector<SOLUTION> my_subopt(char *seq, int delta){
    std::vector<SOLUTION> ret;
    SOLUTION *sol = subopt(seq, NULL, delta, NULL);
    for(int i = 0; sol[i].structure != NULL; i++){
      ret.push_back(sol[i]);
    }
    free(sol);
    /* hopefully the structures will be free'd when the vector is destroyed */
    return ret;
  }
  SOLUTION *my_subopt(char *seq, char *constraint, int delta){
    return subopt(seq, constraint, delta, NULL);
  }

%}

%newobject subopt;

SOLUTION *my_subopt(char *seq, char *constraint, int delta, FILE *fp);
std::vector<SOLUTION> my_subopt(char *seq, int delta, FILE *fp);
std::vector<SOLUTION> my_subopt(char *seq, int delta);
SOLUTION *my_subopt(char *seq, char *constraint, int delta);

%ignore subopt;
%ignore subopt_par;
%ignore subopt_circ;
%ignore zukersubopt;
%ignore zukersubopt_par;

%include  <ViennaRNA/subopt.h>
