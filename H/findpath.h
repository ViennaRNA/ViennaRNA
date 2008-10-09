#ifndef FIND_PATH_H
#define FIND_PATH_H
typedef struct path {
  double en;
  char *s;
} path_t;

extern int find_saddle (char *seq, char *struc1, char *struc2, int max);
extern path_t* get_path(char *seq, char *s1, char* s2, int maxkeep);

#endif
