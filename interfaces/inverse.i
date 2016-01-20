/**********************************************/
/* BEGIN interface for Inverse Folding a.k.a. */
/* sequence design                            */
/**********************************************/

//%subsection "Inverse Folding"

%rename (inverse_fold) my_inverse_fold;
%{
  char *my_inverse_fold(char *start, const char *target, float *cost) {
    char *seq;
    int n;
    n = strlen(target);
    seq = vrna_random_string(n, symbolset);
    if (start)
      strncpy(seq, start, strlen(start));
    *cost = inverse_fold(seq, target);
    if (start)
      /* for backward compatibility modify start */
      strncpy(start, seq, strlen(start));
    return(seq);
  }
%}

%newobject my_inverse_fold;
char * my_inverse_fold(char *start, const char *target, float *OUTPUT);

%rename (inverse_pf_fold) my_inverse_pf_fold;
%{
  char *my_inverse_pf_fold(char *start, const char *target, float *cost) {
    char *seq;
    int n;
    n = strlen(target);
    seq = vrna_random_string(n, symbolset);
    if (start) strncpy(seq, start, n);
    *cost = inverse_pf_fold(seq, target);
    if (start)
      /* for backward compatibility modify start */
      strncpy(start, seq, strlen(start));
    return(seq);
  }
%}

%newobject my_inverse_pf_fold;
char * my_inverse_pf_fold(char *start, const char *target, float *OUTPUT);

%ignore inverse_fold;
%ignore inverse_pf_fold;


%init %{
/* work around segfault when script tries to free symbolset */

symbolset = (char *) vrna_alloc(21);
strcpy(symbolset, "AUGC");

%}

%include  <ViennaRNA/inverse.h>

