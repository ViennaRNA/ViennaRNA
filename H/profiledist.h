
extern float  (*Make_bp_profile(int length))[3];
/* condense pair probability matrix pr into a vector containing probabilities
   for upstream paired, downstream paired and unpaired. This resulting
   probability profile is used as input for profile_edit_distance */

extern float    profile_edit_distance(float (*T1)[3], float (*T2)[3]);
/* align two probability profiles */

extern void     print_bppm(float (*T)[3]);
/* print string representation of probability profile */

extern void     free_profile(float (*T)[3]);
/* free space allocated in Make_bp_profile */
