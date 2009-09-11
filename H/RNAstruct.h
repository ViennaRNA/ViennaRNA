#ifndef __VIENNA_RNA_PACKAGE_RNASTRUCT_H__
#define __VIENNA_RNA_PACKAGE_RNASTRUCT_H__

#define STRUC     1000

char *b2HIT(const char *structure); /* Full   -> HIT    [incl. root] */
char *b2C(const char *structure);   /* Full   -> Coarse [incl. root] */
char *b2Shapiro(const char *structure); /* Full -> weighted Shapiro [i.r.] */
char *add_root(const char *);        /* {Tree} -> ({Tree}R)          */

char  *expand_Shapiro(const char *coarse);
/* add S for stacks to coarse struct */
char  *expand_Full(const char *structure); /* Full   -> FFull         */
char  *unexpand_Full(const char *ffull);   /* FFull  -> Full          */
char  *unweight(const char *wcoarse);   /* remove weights from coarse struct */

void   unexpand_aligned_F(char *align[2]);

void   parse_structure(const char *structure); /* make structure statistics */

extern int    loop_size[STRUC];       /* loop sizes of a structure */
extern int    helix_size[STRUC];      /* helix sizes of a structure */
extern int    loop_degree[STRUC];     /* loop degrees of a structure */
extern int    loops;                  /* n of loops and stacks */
extern int    unpaired, pairs;        /* n of unpaired digits and pairs */

/*
  Example:
  .((..(((...)))..((..)))).   is the bracket or full tree
  becomes expanded:   - expand_Full() -
  ((U)(((U)(U)((((U)(U)(U)P)P)P)(U)(U)(((U)(U)P)P)P)P)(U)R)
  HIT:                - b2HIT() -
  ((U1)((U2)((U3)P3)(U2)((U2)P2)P2)(U1)R) 
  Coarse:             - b2C() -
  ((H)((H)M)R)
  becomes expanded:   - expand_Shapiro() -
  (((((H)S)((H)S)M)S)R)
  weighted Shapiro:   - b2Shapiro() -
  ((((((H3)S3)((H2)S2)M4)S2)E2)R)
*/

#endif
