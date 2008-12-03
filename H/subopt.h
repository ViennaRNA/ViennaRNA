/* subopt.h */
#ifndef SUBOPT_H
#define SUBOPT_H

typedef struct {
  float energy;                            /* energy of structure */
  char *structure;
} SOLUTION;

extern  SOLUTION *subopt (char *seq, char *sequence, int delta, FILE *fp);
extern  SOLUTION *subopt_circ (char *seq, char *sequence, int delta, FILE *fp);
		/* returns list of subopt structures or writes to fp */

extern  int subopt_sorted;                           /* sort output by energy */
#endif
/* End of file */
