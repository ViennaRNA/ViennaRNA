#define NBASES 8
static char Law_and_Order[] = "_ACGUXKI";
static int BP_pair[NBASES][NBASES]=
/* _  A  C  G  U  X  K  I */
{{ 0, 0, 0, 0, 0, 0, 0, 0},
 { 0, 0, 0, 0, 5, 0, 0, 5},
 { 0, 0, 0, 1, 0, 0, 0, 0},
 { 0, 0, 2, 0, 3, 0, 0, 0},
 { 0, 6, 0, 4, 0, 0, 0, 6},
 { 0, 0, 0, 0, 0, 0, 2, 0},
 { 0, 0, 0, 0, 0, 1, 0, 0},
 { 0, 6, 0, 0, 5, 0, 0, 0}};

#define MAXALPHA 20       /* maximal length of alphabet */

static short alias[MAXALPHA+1];
static int pair[MAXALPHA+1][MAXALPHA+1];

#define ENCODE(C) (energy_set>0)?((C)-'A'+1): \
   ((strchr(Law_and_Order, (C))==0)?0:(strchr(Law_and_Order, (C))-Law_and_Order))

extern char *nonstandards;
static void make_pair_matrix(void)
{
   int i,j;
   
   if (energy_set==0) {
      for (i=0; i<5; i++) alias[i] = i;
      alias[5] = 3; /* X <-> G */
      alias[6] = 2; /* K <-> C */
      alias[7] = 0; /* I <-> default base '@' */
      for (i=0; i<NBASES; i++) {
	 for (j=0; j<NBASES; j++) 
	    pair[i][j] = BP_pair[i][j];
      }
      if (noGU) pair[3][4] = pair[4][3] =0;
      if (nonstandards!=NULL) {  /* allow nonstandard bp's */ 
	 for (i=0; i<strlen(nonstandards); i+=2) 
	    pair[ENCODE(nonstandards[i])][ENCODE(nonstandards[i+1])]=7;
      }
   } else {
      for (i=0; i<=MAXALPHA; i++) {
	 for (j=0; j<=MAXALPHA; j++) 
	    pair[i][j] = 0;
      }
      if (energy_set==1) {
	 for (i=1; i<=MAXALPHA;) {
	    alias[i++] = 3;  /* A <-> G */
	    alias[i++] = 2;  /* B <-> C */
	 }
	 for (i=1; i<=MAXALPHA; i++) {
	    pair[i][i+1] = 2;    /* AB <-> GC */
	    i++;
	    pair[i][i-1] = 1;    /* BA <-> CG */
	 }
      }
      else if (energy_set==2) {
	 {
	    alias[i++] = 1;  /* A <-> A*/
	    alias[i++] = 4;  /* B <-> U */
	 }
	 for (i=1; i<=MAXALPHA; i++) {
	    pair[i][i+1] = 5;    /* AB <-> AU */
	    i++;
	    pair[i][i-1] = 6;    /* BA <-> UA */
	 }
      } else nrerror("What energy_set are YOU using??");
   }
}


