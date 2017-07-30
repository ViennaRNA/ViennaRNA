#include <stdio.h>

#include <queue>

#include "barrier_tree.h"

using namespace std;

typedef struct {
  float barrier;
  int i;
  int j;
  bool findpath;
} energy_pair;

struct comparator {
  bool operator()(const energy_pair& x, const energy_pair& y) const {
    if (x.barrier==y.barrier) {
      if (x.findpath != y.findpath) return x.findpath;  // findpath==false - go first!
      if (x.i == y.i) return x.j > y.j;
      else return x.i > y.i;
    }
    return x.barrier > y.barrier;
  }
};

// union-find set array
vector<int> parent;
//unsigned int num_unions = 0;

// ===================== LOCAL FUNCTIONS ====================

// and union-find set functions
int find(int x) {
  if (x != parent[x] && parent[x] != parent[parent[x]])
    parent[x] = find(parent[x]);
  return parent[x];
}

void union_set(int father, int child) {
  int u, v;
  u = find(father);
  v = find(child);
  if (u != v) {
    parent[max(v, u)] = min(v, u);
    //num_unions++;
  }
}

bool joint(int x, int y) {
  return find(x) == find(y);
}

void init_union(int n) {
  parent.resize(n);
  for (int i=0; i<n; i++) parent[i]=i;
}

// find father of node i (uncomputed -> return i)
int findfather(nodeT *nodes, int i) {
  if (nodes[i].father == -1) return i;
  else return findfather(nodes, nodes[i].father);
}

// make barrier tree
int make_tree(int n, float *energy_barr, bool *findpath, nodeT *nodes)
{
  priority_queue<energy_pair, vector<energy_pair>, comparator> saddles;
  for (int i=0; i<n; i++) {
    //if (nodes[i].father!=-1) continue;
    for (int j=i+1; j<n; j++) {
      //if (nodes[j].father!=-1) continue;
      if (energy_barr[i*n+j]<1e8) {
        energy_pair ep;
        ep.barrier = energy_barr[i*n+j];
        ep.i=i;
        ep.j=j;
        ep.findpath=findpath[i*n+j];
        saddles.push(ep);
      }
    }
  }

  // max_height
  float max_height = -1e10;
  //for (int i=0; i<n; i++) if (nodes[i].father!=-1) max_height=max(max_height, nodes[i].saddle_height);

  // compute all except one nodes
  /*while (!saddles.empty()) {
    energy_pair ep = saddles.top();
    saddles.pop();

    // if not already computed
    int fatheri = findfather(nodes, ep.i);
    int fatherj = findfather(nodes, ep.j);
    if (fatheri != fatherj) {
      // merge i and j by their fathers
      if (fatheri>fatherj) swap(fatheri, fatherj);
      nodes[fatherj].saddle_height = ep.barrier;
      add_father(nodes, fatherj, fatheri, 0.5);
      //nodes[fatherj].father = fatheri;
      if (ep.barrier>max_height) max_height = ep.barrier;
    }
  }*/

  float energy = 1e10;
  set<std::pair<int, int> > degen_set;

  init_union(n);
  while (!saddles.empty()) {
    energy_pair ep = saddles.top();
    saddles.pop();

    if (!joint(ep.i, ep.j)) {

      int i=find(ep.i);
      int j=find(ep.j);

      int father = min(i, j);
      int child = max(i, j);

      // degeneracy :/
      if (energy != ep.barrier) {
        degen_set.clear();
      } else {
        for (set<std::pair<int, int> >::iterator it=degen_set.begin(); it!=degen_set.end(); it++) {
          // if last father equals this child - change that father
          if (it->first == child) nodes[it->second].father = father;
        }
      }
      energy = ep.barrier;
      degen_set.insert(make_pair(father, child));

      //fprintf(stderr, "joining: (%2d)%2d (%2d)%2d f%2d ch%2d %6.2f %c\n", ep.i, i, ep.j, j, father, child, ep.barrier, ep.findpath?'#':' ');
      nodes[child].father = father;
      nodes[child].saddle_height = ep.barrier;
      nodes[child].color = (ep.findpath?0.5:0.0);

      if (ep.barrier>max_height) max_height = ep.barrier;

      // finally join them
      union_set(father, child);
    }
  }

  // finish the last one
  nodes[0].saddle_height = max_height + 0.1;
  nodes[0].color = 0.5;

  return 0;
}


void add_father(nodeT *nodes, int child, int father, double color)
{
  //search for old father
  int old_father = nodes[child].father;
  set<int>::iterator it;
  if (old_father !=-1 && (it = nodes[old_father].children.find(child))!=nodes[old_father].children.end()) {
    nodes[old_father].children.erase(it);
  }

  // add new one
  nodes[child].father = father;
  if (color>=0.0) nodes[child].color = color;
  nodes[father].children.insert(child);

  // recompute others
  set<int> tmp = nodes[child].children;
  for (it=tmp.begin(); it!=tmp.end(); it++) {
    if (nodes[*it].saddle_height > nodes[child].saddle_height) {

      //nodes[*it].father = father;
      add_father(nodes, *it, father, -1.0);
    }
  }
}
