/*
  Last changed Time-stamp: <2007-07-10 14:12:49 xtof>
  $Id: options.h,v 1.11 2007/11/03 16:45:58 Kinwalker Exp $
*/
#ifndef _OPTIONS_H_
#define _OPTIONS_H_

#if __cplusplus
#  define BEGIN_C_DECLS extern "C" {
#  define END_C_DECLS   }
#else
#  define BEGIN_C_DECLS
#  define END_C_DECLS
#endif

BEGIN_C_DECLS

typedef struct _Options {

  /*options without argument*/
  char ExeName[256];
  int init_structure;
  int interrupt_trajectory;
  int printfront;
  int testseq;
  int verbose;

  /*options with argument*/
  char barrier_heuristic;
  int dangle;
  char grouping[256];
  int lookahead;

  int maxkeep;
  int noLonelyPairs;
  int transcribed;
  float transcription_rate;
  int windowsize;

  /* int fold_constrained; */
/*   int minLoopSize; */
/*   int minStackSize; */
/*   char Emodel; */

} OptionS;

OptionS *
decodeCML(int argc, char *argv[]);
char*
optionString(void);

END_C_DECLS

#endif /* _OPTIONS_H_ */

/* End of file */
