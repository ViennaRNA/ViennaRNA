#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <string.h>
#include  "utils.h"
#include  "fold_vars.h"
#include  "fold.h"
#include  "part_func.h"
#include  "inverse.h"
#include  "RNAstruct.h"
#include  "treedist.h"
#include  "stringdist.h"
#include  "profiledist.h"

void main()
{
   char *seq1="CGCAGGGAUACCCGCG", *seq2="GCGCCCAUAGGGACGC",
        *struct1,* struct2,* xstruc;
   float e1, e2, tree_dist, string_dist, profile_dist, kT;
   Tree *T1, *T2;
   swString *S1, *S2;
   float *pf1, *pf2;
   FLT_OR_DBL *bppm;
   /* fold at 30C instead of the default 37C */
   temperature = 30.;      /* must be set *before* initializing  */

   /* allocate memory for structure and fold */
   struct1 = (char* ) space(sizeof(char)*(strlen(seq1)+1));
   e1 =  fold(seq1, struct1);

   struct2 = (char* ) space(sizeof(char)*(strlen(seq2)+1));
   e2 =  fold(seq2, struct2);

   free_arrays();     /* free arrays used in fold() */

   /* produce tree and string representations for comparison */
   xstruc = expand_Full(struct1);
   T1 = make_tree(xstruc);
   S1 = Make_swString(xstruc);
   free(xstruc);

   xstruc = expand_Full(struct2);
   T2 = make_tree(xstruc);
   S2 = Make_swString(xstruc);
   free(xstruc);

   /* calculate tree edit distance and aligned structures with gaps */
   edit_backtrack = 1;
   tree_dist = tree_edit_distance(T1, T2);
   free_tree(T1); free_tree(T2);
   unexpand_aligned_F(aligned_line);
   printf("%s\n%s  %3.2f\n", aligned_line[0], aligned_line[1], tree_dist);

   /* same thing using string edit (alignment) distance */
   string_dist = string_edit_distance(S1, S2);
   free(S1); free(S2);
   printf("%s  mfe=%5.2f\n%s  mfe=%5.2f  dist=%3.2f\n",
          aligned_line[0], e1, aligned_line[1], e2, string_dist);

   /* for longer sequences one should also set a scaling factor for
      partition function folding, e.g: */
   kT = (temperature+273.15)*1.98717/1000.;  /* kT in kcal/mol */
   pf_scale = exp(-e1/kT/strlen(seq1));

   /* calculate partition function and base pair probabilities */
   e1 = pf_fold(seq1, struct1);
   /* get the base pair probability matrix for the previous run of pf_fold() */
   bppm = export_bppm();
   pf1 = Make_bp_profile_bppm(bppm, strlen(seq1));

   e2 = pf_fold(seq2, struct2);
   /* get the base pair probability matrix for the previous run of pf_fold() */
   bppm = export_bppm();
   pf2 = Make_bp_profile_bppm(bppm, strlen(seq2));

   free_pf_arrays();  /* free space allocated for pf_fold() */

   profile_dist = profile_edit_distance(pf1, pf2);
   printf("%s  free energy=%5.2f\n%s  free energy=%5.2f  dist=%3.2f\n",
          aligned_line[0], e1, aligned_line[1], e2, profile_dist);

   free_profile(pf1); free_profile(pf2);
}
