/**********************************************/
/* BEGIN interface for Parsing and Comparing  */
/**********************************************/

//%section "Parsing and Comparing Structures"
// from RNAstruct.h

%newobject b2HIT;
%newobject b2C;
%newobject b2Shapiro;
%newobject add_root;
%newobject expand_Shapiro;
%newobject expand_Full;
%newobject unexpand_Full;
%newobject unweight;
char *b2HIT(char *structure);     // Full   -> HIT    [incl. root]
char *b2C(char *structure);       // Full   -> Coarse [incl. root]
char *b2Shapiro(char *structure); // Full -> weighted Shapiro [i.r]
char *add_root(char *);           // {Tree} -> ({Tree}R)
char *expand_Shapiro(char *coarse); // add S for stacks to coarse struct
char *expand_Full(char *structure); // Full   -> FFull
char *unexpand_Full(char *ffull);   // FFull  -> Full
char *unweight(char *wcoarse);   // remove weights from coarse struct
void   unexpand_aligned_F(char *align[2]);
void   parse_structure(char *structure); // make structure statistics
int    loop_size[1000];       // loop sizes of a structure
int    helix_size[1000];      // helix sizes of a structure
int    loop_degree[1000];     // loop degrees of a structure
int    loops;                 // n of loops and stacks
int    unpaired, pairs;       // n of unpaired digits and pairs

%include  <ViennaRNA/treedist.h>
%include  <ViennaRNA/stringdist.h>
%newobject Make_bp_profile;
%include  <ViennaRNA/profiledist.h>
// from dist_vars.h
int   edit_backtrack;  /* set to 1 if you want backtracking */
char *aligned_line[2]; /* containes alignment after backtracking */
int  cost_matrix;      /* 0 usual costs (default), 1 Shapiro's costs */

// this doesnt work currently
%inline %{
void *deref_any(void **ptr, int index) {
   /* dereference arbitray pointer */
   return (void *) ptr[index];
}
%}

%{
char *get_aligned_line(int i) {
  i = i % 2;
  return aligned_line[i];
}
%}

char *get_aligned_line(int);

