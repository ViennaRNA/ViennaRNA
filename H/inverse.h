/* prototypes from inverse.c */
extern char *symbolset;    /* alphabet default is "AUGC" */
extern float inverse_fold(char *start, const char *target);  
/* find sequences with predefined structure.
   the found sequence is written to start,
   return value is
      energy_of_struct(start, target) - fold(start, structure),
   i.e. 0. if search was successful; */
   
extern float inverse_pf_fold(char *start, const char *target);
/*  inverse folding maximising the frequency of target in the
    ensemble of structures, final sequence is written to start, returns 
       energy_of_struct(start, target) - part_func(start, structure)
*/
extern float final_cost; /* when to stop inverse_pf_fold() */
extern int give_up;   /* default 0: try to minimize structure distance even if 
			 no exact solution can be found */

