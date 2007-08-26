#ifndef DUPLEX_H
#define DUPLEX_H

typedef struct {
  int i;
  int j;
  char *structure;
  float energy;
} duplexT;

extern duplexT duplexfold(const char *s1, const char *s2);
extern duplexT *duplex_subopt(const char *s1, const char *s2, int delta, int w);
extern duplexT aliduplexfold(const char *s1[], const char *s2[]);
extern duplexT *aliduplex_subopt(const char *s1[], const char *s2[], int delta, int w);

#endif
