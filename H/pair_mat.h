#include <ctype.h>
#define NBASES 8
/*@notnull@*/

static const char Law_and_Order[] = "_ACGUTXKI";
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
/* rtype[pair[i][j]]:=pair[j][i] */
static int rtype[8] = {0, 2, 1, 4, 3, 6, 5, 7};

#ifdef _OPENMP
#pragma omp threadprivate(Law_and_Order, BP_pair, alias, pair, rtype)
#endif

/* for backward compatibility */
#define ENCODE(c) encode_char(c)

static int encode_char(char c) {
  /* return numerical representation of base used e.g. in pair[][] */
  int code;
  if (energy_set>0) code = (int) (c-'A')+1;
  else {
    char *pos;
    pos = strchr(Law_and_Order, c);
    if (pos==NULL) code=0;
    else code = (int) (pos-Law_and_Order);
    if (code>5) code = 0;
    if (code>4) code--; /* make T and U equivalent */
  }
  return code;
}

/*@+boolint +charint@*/
/*@null@*/
extern char *nonstandards;
extern void   nrerror(const char message[]);
static void make_pair_matrix(void)
{
   int i,j;

   if (energy_set==0) {
      for (i=0; i<5; i++) alias[i] = (short) i;
      alias[5] = 3; /* X <-> G */
      alias[6] = 2; /* K <-> C */
      alias[7] = 0; /* I <-> default base '@' */
      for (i=0; i<NBASES; i++) {
          for (j=0; j<NBASES; j++)
            pair[i][j] = BP_pair[i][j];
      }
      if (noGU) pair[3][4] = pair[4][3] =0;
      if (nonstandards!=NULL) {  /* allow nonstandard bp's */
         for (i=0; i<(int)strlen(nonstandards); i+=2)
            pair[encode_char(nonstandards[i])]
              [encode_char(nonstandards[i+1])]=7;
      }
      for (i=0; i<NBASES; i++) {
          for (j=0; j<NBASES; j++)
           rtype[pair[i][j]] = pair[j][i];
      }
   } else {
      for (i=0; i<=MAXALPHA; i++) {
         for (j=0; j<=MAXALPHA; j++)
            pair[i][j] = 0;
      }
      if (energy_set==1) {
         for (i=1; i<MAXALPHA;) {
            alias[i++] = 3;  /* A <-> G */
            alias[i++] = 2;  /* B <-> C */
         }
         for (i=1; i<MAXALPHA; i++) {
            pair[i][i+1] = 2;    /* AB <-> GC */
            i++;
            pair[i][i-1] = 1;    /* BA <-> CG */
         }
      }
      else if (energy_set==2) {
        for (i=1; i<MAXALPHA;) {
            alias[i++] = 1;  /* A <-> A*/
            alias[i++] = 4;  /* B <-> U */
         }
         for (i=1; i<MAXALPHA; i++) {
            pair[i][i+1] = 5;    /* AB <-> AU */
            i++;
            pair[i][i-1] = 6;    /* BA <-> UA */
         }
      }
      else if (energy_set==3) {
        for (i=1; i<MAXALPHA-2; ) {
          alias[i++] = 3;  /* A <-> G */
          alias[i++] = 2;  /* B <-> C */
          alias[i++] = 1;  /* C <-> A */
          alias[i++] = 4;  /* D <-> U */
        }
        for (i=1; i<MAXALPHA-2; i++) {
          pair[i][i+1] = 2;    /* AB <-> GC */
          i++;
          pair[i][i-1] = 1;    /* BA <-> CG */
          i++;
          pair[i][i+1] = 5;    /* CD <-> AU */
          i++;
          pair[i][i-1] = 6;    /* DC <-> UA */
        }
      }
      else nrerror("What energy_set are YOU using??");
      for (i=0; i<=MAXALPHA; i++) {
        for (j=0; j<=MAXALPHA; j++)
          rtype[pair[i][j]] = pair[j][i];
      }
   }
}

static short *encode_sequence(const char *sequence, short how){
  unsigned int i,l = (unsigned int)strlen(sequence);
  short         *S = (short *) space(sizeof(short)*(l+2));

  switch(how){
    /* standard encoding as always used for S */
    case 0:   for(i=1; i<=l; i++) /* make numerical encoding of sequence */
                S[i]= (short) encode_char(toupper(sequence[i-1]));
              S[l+1] = S[1];
              S[0] = (short) l;
              break;
    /* encoding for mismatches of nostandard bases (normally used for S1) */
    case 1:   for(i=1; i<=l; i++)
                S[i] = alias[(short) encode_char(toupper(sequence[i-1]))];
              S[l+1] = S[1];
              S[0] = S[l];
              break;
  }

  return S;
}
