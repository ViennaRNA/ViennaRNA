#ifndef PKPLEX_H
#define PKPLEX_H

/* remove or comment the typedefs for SOLUTION2 and duplexT if you get "conflicting types" error messages during compilation. This should only happen if your Vienna RNA package version is newer than the one available at the release of PKplex. */

/*typedef struct {
  float energy;                            
  char *structure;
  int i;
  int j;
} SOLUTION2;
typedef struct {
  int i;
  int j;
  char *structure;
  float energy;
} duplexT;*/

typedef struct dupVar{
  int i;
  int j;
  int end;
  char *structure;
  double energy;
  int offset;
  double dG1;
  double dG2;
  double ddG;
  int tb;
  int te;
  int qb;
  int qe;
  //} duplexT;
} dupVar;

extern int subopt_sorted;

extern dupVar** PKLduplexfold_XS(const char*s1, int **access_s1, const int threshold, const int alignment_length, const int delta);//, const int target_dead, const int query_dead);

extern int      arraySize(duplexT** array);
extern void     freeDuplexT(duplexT** array);

#endif
