extern char symbolset[];    /* alphabet default is "AUGC" */
/* the found sequence is given back in start,
   return value for inverse_fold is
          energy_of_struct(start, target) - fold(start, structure),
   that is 0. if search was successful;
   for inverse_pf_fold
          energy_of_struct(start, target) - part_func(start, structure)
   that is the frequency of target in the ensemble of strucures is
          p = exp(-inverse_pf_fold/kT);
*/
extern float inverse_fold(char *start, char *target);  
extern float inverse_pf_fold(char *start, char *target);
extern float final_cost; /* when to stop inverse_pf_fold() */
extern int give_up;   /* default 0: try to minimize structure distance even if 
			 no exact solution can be found */

