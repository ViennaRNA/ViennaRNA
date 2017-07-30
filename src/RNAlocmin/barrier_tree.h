#include "treeplot.h"


// union find set for LM when trying to recompute barrier tree
void union_set(int father, int child);
void init_union(int n);
int find(int x);

// make barrier tree
int make_tree(int n, float *energy_bar, bool *findpath, nodeT *nodes);

// recompute single father change
void add_father(nodeT *nodes, int child, int father, double color);
