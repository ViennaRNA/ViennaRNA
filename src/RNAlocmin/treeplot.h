#ifndef __TREEPLOT_H
#define __TREEPLOT_H

#include <set>

using namespace std;

/* treeplot.h */
typedef struct {
  float height;         /* height (energy, time, whatever) of this leaf    */
  float saddle_height;  /* height of internal node that connects this leaf */
  int father;           /* node with which it connects                     */
  float color;            /* color of the connection 0=black 1=white         */
  char *label;          /* label string, if NULL use index+1               */
  set<int> children;        /* set of children */
} nodeT;

// generate barrier tree (from n nodes to filename) (stolen from barriers)
void PS_tree_plot(nodeT *nodes, int n, char *filename);

#endif
/* End of file */
