//#include <vector>
//#include <iostream>
//using namespace std;
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>
#include <random>
#include <stack>
#include <iomanip>
#include <cstdint>
#include <vector>
#include <iostream>

#ifdef _WIN32
#ifdef __MINGW32__
#include <unistd.h>
#else
#include "ViennaRNA/intern/unistd_win.h"
#endif
#else
#include <unistd.h>
#endif

using namespace std;


#ifdef __cplusplus
extern "C" {
#endif
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "string.h"
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#ifdef __cplusplus
}
#endif

#include "linearbpseeds.h" //#include "ViennaRNA/inverse/linearbpseeds.h"


static random_device                rd;
static mt19937                      gen(rd());
static uniform_real_distribution<>  distrib(0.0, 1.0);

//There are 8 possibilities per node assignment mapped to 8 different int: GC <-> 0, CG <-> 1, AU <-> 2, UA <-> 3, ?? <-> 4, A <-> 5, C <-> 6, RR (i.e root) <-> 7
//get_assignments(_biseparated) is optimized only for multiloops with a degree below 64 for now.
//Warning: still some risks of overflow with how the states are defined for now both in cache and cacheassignment
static long long                    max_degree_optimized_multiloop = 64;

static constexpr int                delta[8]  = { 1, -1, 0, 0, 0, 0, 0, 0 };
static constexpr int                deltaC[8] = { 0, 0, -1, 1, 0, 0, 0, 0 };

static unordered_map <int, vector<int> >
                                    children_colors_from_parent =
{
  //Possible children assignment given the parent assignment
  { 2, { 2, 0, 1 } },
  { 3, { 3, 0, 1 } },
  { 0, { 3, 0, 2 } },
  { 1, { 3, 1, 2 } },
  { 7, { 0, 3, 1, 2} }
};


struct Tree {
  pair <int, int> parent;
  bool            unproper; //A boolean to know if the loop is a forbidden motif
  vector<Tree *>  children;
  uint64_t        leafbits; //A signature of the children: if they are leaves or internal nodes
  int             nb_children;
  bool            is_leaf;
};


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

static long long
make_leafbits(const vector<Tree *>& l);


static pair<string, string>
unwrap(int c);


static vector <vector<int> >
get_assignments(const vector<Tree *>& l,
                int                   c,
                size_t                i);


static vector <vector<int> >
get_assignments_biseparated(const vector<Tree *>& l,
                            int                   c,
                            int                   color,
                            size_t                i);


static long double
num_design(const Tree                                         *v,
           int                                                c,
           int                                                current_level_mod,
           int                                                m,
           const unordered_set<int>&                          leaves_levels_mod,
           unordered_map<long long, long double>&             cache,
           unordered_map<long long, vector <vector<int> > >&  cacheassignment,
           bool                                               flexible,
           long double                                        GCweight,
           int                                                n);


static long double
num_design_biseparated(const Tree                                         *v,
                       int                                                c,
                       int                                                current_level_modA,
                       int                                                mA,
                       const unordered_set<int>&                          leaves_levels_modA,
                       int                                                current_level_modC,
                       int                                                mC,
                       const unordered_set<int>&                          leaves_levels_modC,
                       unordered_map<long long, long double>&             cache,
                       unordered_map<long long, vector <vector<int> > >&  cacheassignment,
                       bool                                               flexible,
                       long double                                        GCweight,
                       int                                                n);


static string
stochastic_backtrack(const Tree                                         *v,
                     int                                                c,
                     int                                                current_level_mod,
                     int                                                m,
                     const unordered_set<int>&                          leaves_levels_mod,
                     unordered_map<long long, long double>&             cache,
                     unordered_map<long long, vector <vector<int> > >&  cacheassignment,
                     bool                                               flexible,
                     long double                                        GCweight,
                     int                                                n);


static string
stochastic_backtrack_biseparated(const Tree *v,
                                 int        c,
                                 int
                                 current_level_modA,
                                 int mA,
                                 const unordered_set<int>&
                                 leaves_levels_modA,
                                 int
                                 current_level_modC,
                                 int mC,
                                 const unordered_set<int>&
                                 leaves_levels_modC,
                                 unordered_map<long long, long double>&             cache,
                                 unordered_map<long long, vector <vector<int> > >&  cacheassignment,
                                 bool                                               flexible,
                                 long double                                        GCweight,
                                 int                                                n);


static vector<int>
ssparse(const string& seq);


static bool
checker(const vector<Tree *>& v,
        int                   lim1,
        int                   lim2);


static vector<Tree *>
dbn_to_tree_list(const vector<int>& dbn,
                 int                i,
                 int                j);


static Tree *
dbn_to_tree(const vector<int>& dbn);


static vector<vector<int> >
part(int i);


static void
tree_cleanup(Tree *v);


static unordered_set<string>
sample(size_t                                             magic_threshold,
       long double                                        Zfull,
       const Tree                                         *tree,
       int                                                m,
       vector<unordered_set<int> >                        tuple_leaf_level_list,
       vector<long double>                                Zlist,
       vector<unordered_map<long long, long double> >     cachelist,
       unordered_map<long long, vector <vector<int> > >&  cacheassignment,
       bool                                               flex,
       long double                                        GCweight,
       int                                                n);


static unordered_set<string>
sample_biseparated(size_t                                             magic_threshold,
                   long double                                        Zfull,
                   const Tree                                         *tree,
                   int                                                mA,
                   int                                                mC,
                   vector<unordered_set<int> >                        tuple_leaf_levelA_list,
                   vector<unordered_set<int> >                        tuple_leaf_levelC_list,
                   vector<long double>                                Zlist,
                   vector<unordered_map<long long, long double> >     cachelist,
                   unordered_map<long long, vector <vector<int> > >&  cacheassignment,
                   bool                                               flex,
                   long double                                        GCweight,
                   int                                                n);


vector<string>
generate_seeds(string                       structure,
               unsigned int                 target_nb_sols,
               vrna_inverse_linearbp_opt_t  *p);


/*
 #########################################
 # BEGIN OF PUBLIC FUNCTION DEFINITIONS  #
 #      (all available in RNAlib)        #
 #########################################
 */
PUBLIC char **
vrna_inverse_seed_linearbp(const char                   *structure,
                           size_t                       num_samples,
                           vrna_inverse_linearbp_opt_t  *options)
{
  char                        **samples   = NULL;
  size_t                      ns          = 0;
  vrna_inverse_linearbp_opt_t default_opt = {
    (unsigned int)0,
    (unsigned int)3,
    (double)0.0
  };

  if (structure) {
    if (options == NULL)
      options = &default_opt;

    vector<string> res = generate_seeds(string(structure),
                                        num_samples,
                                        options);

    if (res.size() > 0) {
      ns      = res.size();
      samples = (char **)vrna_alloc(sizeof(char *) * (ns + 1));

      for (size_t i = 0; i < ns; i++)
        samples[i] = strdup(res[i].c_str());

      samples[ns] = NULL; /* end of list indicator */
    }
  }

  return samples;
}


/*
 #########################################
 # BEGIN OF PRIVATE FUNCTION DEFINITIONS #
 #          (internal use only)          #
 #########################################
 */
static long long
make_leafbits(const vector<Tree *>& l)
{
  /*
   *  Input:
   * l, a list of secondary structure tree
   *  Output:
   * A boolean, 1 if the node v->parent is a leaf or not
   */
  long long leafbits = 0;

  for (size_t i = 0; i < l.size(); i++)
    if ((l[i])->is_leaf)
      leafbits |= (uint64_t(1) << i);

  return leafbits;
}


static pair<string, string>
unwrap(int c)
{
  /*
   *  Input:
   * c, a base pair as an int
   *  Output:
   * A couple of characters corresponding to the base pair
   */
  if (c == 0) {
    return { "G", "C" };
  }

  if (c == 1) {
    return { "C", "G" };
  }

  if (c == 2) {
    return { "A", "U" };
  }

  if (c == 3) {
    return { "U", "A" };
  }

  throw invalid_argument("Should not reach here: not a valid base pair");
}


static vector <vector<int> >
get_assignments(const vector<Tree *>& l,
                int                   c,
                size_t                i)
{
  /*
   *  Input:
   * l, an ordered list of the trees for the leaves and base pair children
   * c, assignment of the parent base pair (or information that the parent is the root)
   * i, the index of the node considered to assign
   *  Output:
   * The complete list of the ordered list of assignments for the leaves and base pair children with A at unpaired positions
   */
  if (i >= l.size()) {
    return { {} };
  } else {
    Tree *v = l[i];
    if (v->is_leaf) {
      if ((c != 2) && (c != 3)) {
        const auto&           temp = get_assignments(l, c, i + 1);
        //vector <vector<int>> temp = get_assignments(l,c,i+1);
        vector <vector<int> > resu = {};
        for (size_t i_temp = 0; i_temp < temp.size(); i_temp++) {
          vector<int> temploc     = temp[i_temp];
          uint8_t     used_color  = 0;
          for (size_t j = 0; j < temploc.size(); j++)
            used_color |= (uint8_t(1) << temploc[j]);
          if ((!(used_color & (uint8_t(1) << 2))) && (!(used_color & (uint8_t(1) << 3)))) {
            temploc.insert(temploc.begin(), 5);
            resu.emplace_back(temploc);
          }
        }
        return resu;
      } else {
        return {};
      }
    } else {
      vector <vector<int> > resu  = {};
      const auto&           temp  = get_assignments(l, c, i + 1);
      //vector <vector<int>> temp = get_assignments(l,c,i+1);
      for (size_t cv_i = 0; cv_i < children_colors_from_parent[c].size(); cv_i++) {
        int cv = children_colors_from_parent[c][cv_i];
        for (size_t i_a = 0; i_a < temp.size(); i_a++) {
          vector<int> a           = temp[i_a];
          uint8_t     used_color  = 0;
          for (size_t j = 0; j < a.size(); j++)
            used_color |= (uint8_t(1) << a[j]);
          if (!(used_color & (uint8_t(1) << cv))) {
            if (((cv != 2) && (cv != 3)) || (!(used_color & (uint8_t(1) << 5)))) {
              a.insert(a.begin(), cv);
              resu.emplace_back(a);
            }
          }
        }
      }
      return resu;
    }
  }
}


static vector <vector<int> >
get_assignments_biseparated(const vector<Tree *>& l,
                            int                   c,
                            int                   color,
                            size_t                i)
{
  /*
   *  Input:
   * l, an ordered list of the trees for the leaves and BP children
   * c, assignment of the parent base pair (or information that the parent is the root)
   * color, 0 <-> leaves not determined, 1 <-> leaves are all As, 2 <-> leaves are all Cs
   * i, the index of the considered node to assign
   *  Output:
   * The complete list of the ordered list of assignments for the leaf and base pair children with both A and C at unpaired positions
   */
  if (i >= l.size()) {
    return { {} };
  } else {
    Tree                  *v    = l[i];
    vector <vector<int> > resu  = {};
    if (v->is_leaf) {
      if ((c != 0) && (c != 1) && (color != 1)) {
        const auto& temp = get_assignments_biseparated(l, c, 2, i + 1);
        //vector <vector<int>> temp = get_assignments_biseparated(l,c,2, i+1);
        for (size_t i_temp = 0; i_temp < temp.size(); i_temp++) {
          vector<int> temploc     = temp[i_temp];
          uint8_t     used_color  = 0;
          for (size_t j = 0; j < temploc.size(); j++)
            used_color |= (uint8_t(1) << temploc[j]);
          if ((!(used_color & (uint8_t(1) << 0))) && (!(used_color & (uint8_t(1) << 1)))) {
            temploc.insert(temploc.begin(), 6);
            resu.emplace_back(temploc);
          }
        }
      }

      if ((c != 2) && (c != 3) && (color != 2)) {
        const auto& temp = get_assignments_biseparated(l, c, 1, i + 1);
        //vector <vector<int>> temp = get_assignments_biseparated(l,c,1, i+1);
        for (size_t i_temp = 0; i_temp < temp.size(); i_temp++) {
          vector<int> temploc     = temp[i_temp];
          uint8_t     used_color  = 0;
          for (size_t j = 0; j < temploc.size(); j++)
            used_color |= (uint8_t(1) << temploc[j]);
          if ((!(used_color & (uint8_t(1) << 2))) && (!(used_color & (uint8_t(1) << 3)))) {
            temploc.insert(temploc.begin(), 5);
            resu.emplace_back(temploc);
          }
        }
      }
    } else {
      vector <vector<int> > temp = get_assignments_biseparated(l, c, color, i + 1);
      for (size_t cv_i = 0; cv_i < children_colors_from_parent[c].size(); cv_i++) {
        int cv = children_colors_from_parent[c][cv_i];
        for (size_t i_a = 0; i_a < temp.size(); i_a++) {
          vector<int> a           = temp[i_a];
          uint8_t     used_color  = 0;
          for (size_t j = 0; j < a.size(); j++)
            used_color |= (uint8_t(1) << a[j]);
          if (!(used_color & (uint8_t(1) << cv))) {
            bool boo = 1;
            if (((cv == 2) || (cv == 3)) && (used_color & (uint8_t(1) << 5))) {
              boo = 0;
            }

            if (((cv == 0) || (cv == 1)) && (used_color & (uint8_t(1) << 6))) {
              boo = 0;
            }

            if (boo == 1) {
              a.insert(a.begin(), cv);
              resu.emplace_back(a);
            }
          }
        }
      }
    }

    return resu;
  }
}


static long double
num_design(const Tree                                         *v,
           int                                                c,
           int                                                current_level_mod,
           int                                                m,
           const unordered_set<int>&                          leaves_levels_mod,
           unordered_map<long long, long double>&             cache,
           unordered_map<long long, vector <vector<int> > >&  cacheassignment,
           bool                                               flexible,
           long double                                        GCweight,
           int                                                n)
{
  /*
   *  Input:
   * v, a secondary structure tree with a boolean stored in each node for loop assignment 0 <-> proper, 1 <-> unproper
   * c, the assignment for node v[0], 8 possibilities attributes to 8 different int: GC <-> 0, CG <-> 1, AU <-> 2, UA <-> 3, ?? <-> 4, A <-> 5, C <-> 6, RR <-> 7
   * current_level_mod, the current level modulo m
   * m, the modulus considered for m-separability
   * leaves_levels_mod, the set of levels specific to the leaves
   * cache, the partial solutions for subtrees
   * cacheassignment, the list of local assignments for loops already encountered
   * flexible, a boolean put to 1 if unproper should be the default behavior for every loop
   * GCweight, the weight to attribute to GC base pairs, the default choice should be 0.0
   * n, size of the sequence
   *  Output:
   * The number of separated designs for the current tree under these parameters
   */
  long long state = v->parent.first * (n + 1) * 8 * m + v->parent.second * 8 * m + c * m +
                    current_level_mod;                                                                   //root is assigned (0, n), usual max for pair is (n - 1, n - 1)

  if (!cache.count(state)) {
    if (c == 4) {
      if (v->is_leaf) {
        cache[state] = num_design(v,
                                  5,
                                  current_level_mod,
                                  m,
                                  leaves_levels_mod,
                                  cache,
                                  cacheassignment,
                                  flexible,
                                  GCweight,
                                  n);
      } else {
        cache[state] = 0.0;
        for (size_t couple = 0; couple < 4; couple++)
          cache[state] += num_design(v,
                                     couple,
                                     current_level_mod,
                                     m,
                                     leaves_levels_mod,
                                     cache,
                                     cacheassignment,
                                     flexible,
                                     GCweight,
                                     n);
      }

      return cache[state];
    }

    if (v->is_leaf) {
      if (leaves_levels_mod.count(current_level_mod)) {
        cache[state] = 1.0;
      } else {
        cache[state] = 0.0;
      }
    } else if ((c == 2 || c == 3) && (leaves_levels_mod.count(current_level_mod))) {
      cache[state] = 0.0;
    } else {
      int         next_level_mod = (m + (current_level_mod + delta[c]) % m) % m;
      long double acc;
      if ((!(v->unproper)) && (!flexible)) {
        acc = 0.0;
        long long k = (v->leafbits) * max_degree_optimized_multiloop * 8 + (v->nb_children) * 8 +
                      c;
        if (!cacheassignment.count(k)) {
          cacheassignment[k] = get_assignments(v->children, c, 0);
        }

        const auto& assignments = cacheassignment[k];
        //vector <vector<int>> assignments = cacheassignment[k];
        if (k < 0) {
          cacheassignment.erase(k);
        }

        for (size_t assignment_index = 0; assignment_index < assignments.size();
             assignment_index++) {
          long double prod = 1.0;
          for (size_t i = 0; i < (v->children).size(); i++)
            prod = prod * num_design((v->children)[i],
                                     assignments[assignment_index][i],
                                     next_level_mod,
                                     m,
                                     leaves_levels_mod,
                                     cache,
                                     cacheassignment,
                                     flexible,
                                     GCweight,
                                     n);
          acc += prod;
        }
      } else {
        acc = 1.0;
        for (size_t i = 0; i < (v->children).size(); i++)
          acc *= num_design((v->children)[i],
                            4,
                            next_level_mod,
                            m,
                            leaves_levels_mod,
                            cache,
                            cacheassignment,
                            flexible,
                            GCweight,
                            n);
      }

      if ((GCweight != 0.0) && (c == 0 || c == 1)) {
        acc *= exp(GCweight);
      }

      cache[state] = acc;
    }
  }

  return cache[state];
}


static long double
num_design_biseparated(const Tree                                         *v,
                       int                                                c,
                       int                                                current_level_modA,
                       int                                                mA,
                       const unordered_set<int>&                          leaves_levels_modA,
                       int                                                current_level_modC,
                       int                                                mC,
                       const unordered_set<int>&                          leaves_levels_modC,
                       unordered_map<long long, long double>&             cache,
                       unordered_map<long long, vector <vector<int> > >&  cacheassignment,
                       bool                                               flexible,
                       long double                                        GCweight,
                       int                                                n)
{
  /*
   * Input:
   * v, a secondary structure tree with a boolean stored in each node for loop assignment 0 <-> proper, 1 <-> unproper
   * c, the assignment for node v[0], 8 possibilities attributes to 8 different int: GC <-> 0, CG <-> 1, AU <-> 2, UA <-> 3, ?? <-> 4, A <-> 5, C <-> 6, RR <-> 7
   * current_level_modA (resp C), the current level modulus mA (resp. mC)
   * mA and mC, the modulo considered for (mA, mC)-separability
   * leaves_levels_modA (resp. C), the set of levels specific to the leaves assigned A (resp. C)
   * cache, the partial solutions for subtrees
   * cacheassignment, the list of local assignments for loops already encountered
   * flexible, a boolean put to 1 if unproper should be the default behavior for every loop
   * GCweight, the weight to attribute to GC base pairs, the default choice should be 0.0
   * n, size of the sequence
   * Output:
   * The number of separated designs for the current tree under these parameters
   */
  long long state = v->parent.first * (n + 1) * 8 * mA * mC + v->parent.second * 8 * mA * mC + c *
                    mA * mC + current_level_modA * mC + current_level_modC;

  if (!cache.count(state)) {
    if (c == 4) {
      if (v->is_leaf) {
        cache[state] = num_design_biseparated(v,
                                              5,
                                              current_level_modA,
                                              mA,
                                              leaves_levels_modA,
                                              current_level_modC,
                                              mC,
                                              leaves_levels_modC,
                                              cache,
                                              cacheassignment,
                                              flexible,
                                              GCweight,
                                              n) + num_design_biseparated(v,
                                                                          6,
                                                                          current_level_modA,
                                                                          mA,
                                                                          leaves_levels_modA,
                                                                          current_level_modC,
                                                                          mC,
                                                                          leaves_levels_modC,
                                                                          cache,
                                                                          cacheassignment,
                                                                          flexible,
                                                                          GCweight,
                                                                          n);
      } else {
        cache[state] = 0.0;
        for (size_t couple = 0; couple < 4; couple++)
          cache[state] += num_design_biseparated(v,
                                                 couple,
                                                 current_level_modA,
                                                 mA,
                                                 leaves_levels_modA,
                                                 current_level_modC,
                                                 mC,
                                                 leaves_levels_modC,
                                                 cache,
                                                 cacheassignment,
                                                 flexible,
                                                 GCweight,
                                                 n);
      }

      return cache[state];
    }

    if (v->is_leaf) {
      if ((c == 5) && (leaves_levels_modA.count(current_level_modA))) {
        cache[state] = 1.0;
      } else if ((c == 6) && (leaves_levels_modC.count(current_level_modC))) {
        cache[state] = 1.0;
      } else {
        cache[state] = 0.0;
      }
    } else if ((c == 2 || c == 3) && (leaves_levels_modA.count(current_level_modA))) {
      cache[state] = 0.0;
    } else if ((c == 0 || c == 1) && (leaves_levels_modC.count(current_level_modC))) {
      cache[state] = 0.0;
    } else {
      int         next_level_modA = (mA + (current_level_modA + delta[c]) % mA) % mA;
      int         next_level_modC = (mC + (current_level_modC + deltaC[c]) % mC) % mC;
      long double acc;
      if ((!(v->unproper)) && (!flexible)) {
        acc = 0.0;
        long long k = (v->leafbits) * max_degree_optimized_multiloop * 8 + (v->nb_children) * 8 +
                      c;
        if (!cacheassignment.count(k)) {
          cacheassignment[k] = get_assignments_biseparated(v->children, c, 0, 0);
        }

        const auto& assignments = cacheassignment[k];
        //vector <vector<int>> assignments = cacheassignment[k];
        if (k < 0) {
          cacheassignment.erase(k);
        }

        for (size_t assignment_index = 0; assignment_index < assignments.size();
             assignment_index++) {
          long double prod = 1.0;
          for (size_t i = 0; i < (v->children).size(); i++)
            prod = prod * num_design_biseparated((v->children)[i],
                                                 assignments[assignment_index][i],
                                                 next_level_modA,
                                                 mA,
                                                 leaves_levels_modA,
                                                 next_level_modC,
                                                 mC,
                                                 leaves_levels_modC,
                                                 cache,
                                                 cacheassignment,
                                                 flexible,
                                                 GCweight,
                                                 n);
          acc += prod;
        }
      } else {
        acc = 1.0;
        for (size_t i = 0; i < (v->children).size(); i++)
          acc *= num_design_biseparated((v->children)[i],
                                        4,
                                        next_level_modA,
                                        mA,
                                        leaves_levels_modA,
                                        next_level_modC,
                                        mC,
                                        leaves_levels_modC,
                                        cache,
                                        cacheassignment,
                                        flexible,
                                        GCweight,
                                        n);
      }

      if ((GCweight != 0.0) && (c == 0 || c == 1)) {
        acc *= exp(GCweight);
      }

      cache[state] = acc;
    }
  }

  return cache[state];
}


static string
stochastic_backtrack(const Tree                                         *v,
                     int                                                c,
                     int                                                current_level_mod,
                     int                                                m,
                     const unordered_set<int>&                          leaves_levels_mod,
                     unordered_map<long long, long double>&             cache,
                     unordered_map<long long, vector <vector<int> > >&  cacheassignment,
                     bool                                               flexible,
                     long double                                        GCweight,
                     int                                                n)
{
  /*
   *  Input:
   * v, a secondary structure tree with a boolean stored in each node for loop assignment 0 <-> proper, 1 <-> unproper
   * c, the assignment for node v[0], 8 possibilities attributes to 8 different int: GC <-> 0, CG <-> 1, AU <-> 2, UA <-> 3, ?? <-> 4, A <-> 5, C <-> 6, RR <-> 7
   * current_level_mod, the current level modulo m
   * m, the modulus considered for m-separability
   * leaves_levels_mod, the set of levels specific to the leaves
   * cache, the partial solutions for subtrees
   * cacheassignment, the list of local assignments for loops already encountered
   * flexible, a boolean put to 1 if unproper should be the default behavior for every loop
   * GCweight, the weight to attribute to GC base pairs, the default choice should be 0.0
   * n, size of the sequence
   *  Output:
   * A design uniformly sampled for the current tree under these parameters
   */
  long long state = v->parent.first * (n + 1) * 8 * m + v->parent.second * 8 * m + c * m +
                    current_level_mod;

  if (!cache.count(state)) {
    throw runtime_error("cache miscomputed");
  }

  if (v->is_leaf) {
    return "A";      //factorized on top, no need of specific cases
  }

  if (c == 4) {
    long double r = distrib(gen) * cache[state];
    for (size_t couple = 0; couple < 4; couple++) {
      r -= num_design(v,
                      couple,
                      current_level_mod,
                      m,
                      leaves_levels_mod,
                      cache,
                      cacheassignment,
                      flexible,
                      GCweight,
                      n);
      if (r < 0.0) {
        return stochastic_backtrack(v,
                                    couple,
                                    current_level_mod,
                                    m,
                                    leaves_levels_mod,
                                    cache,
                                    cacheassignment,
                                    flexible,
                                    GCweight,
                                    n);
      }
    }
  } else {
    int next_level_mod = (m + (current_level_mod + delta[c]) % m) % m;
    if ((!(v->unproper)) && (!flexible)) {
      long double r = distrib(gen) * cache[state];
      long long   k = (v->leafbits) * max_degree_optimized_multiloop * 8 + (v->nb_children) * 8 + c;
      if (!cacheassignment.count(k)) {
        throw runtime_error("cache miscomputed");
      }

      const auto& assignments = cacheassignment[k];
      //vector <vector<int>> assignments = cacheassignment[k];
      for (size_t assignment_index = 0; assignment_index < assignments.size(); assignment_index++) {
        long double prod      = 1.0;
        long double GCfactor  = exp(GCweight);
        if ((GCweight != 0.0) && (c == 0 || c == 1)) {
          prod = GCfactor;
        }

        for (size_t i = 0; i < (v->children).size(); i++)
          prod *= num_design((v->children)[i],
                             assignments[assignment_index][i],
                             next_level_mod,
                             m,
                             leaves_levels_mod,
                             cache,
                             cacheassignment,
                             flexible,
                             GCweight,
                             n);
        r -= prod;
        if (r < 0.0) {
          string res = "";
          for (size_t i = 0; i < (v->children).size(); i++) {
            string loc = stochastic_backtrack((v->children)[i],
                                              assignments[assignment_index][i],
                                              next_level_mod,
                                              m,
                                              leaves_levels_mod,
                                              cache,
                                              cacheassignment,
                                              flexible,
                                              GCweight,
                                              n);
            res += loc;
          }
          if (c == 7) {
            return res;
          } else {
            pair <string, string> p = unwrap(c);
            return p.first + res + p.second;
          }
        }
      }
    } else {
      string res = "";
      for (size_t i = 0; i < (v->children).size(); i++) {
        string loc = stochastic_backtrack((v->children)[i],
                                          4,
                                          next_level_mod,
                                          m,
                                          leaves_levels_mod,
                                          cache,
                                          cacheassignment,
                                          flexible,
                                          GCweight,
                                          n);
        res += loc;
      }
      if (c == 7) {
        return res;
      } else {
        pair <string, string> p = unwrap(c);
        return p.first + res + p.second;
      }
    }
  }

  throw invalid_argument("Should not reach here, stochastic backtrack failed");
}


static string
stochastic_backtrack_biseparated(const Tree *v,
                                 int        c,
                                 int
                                 current_level_modA,
                                 int mA,
                                 const unordered_set<int>&
                                 leaves_levels_modA,
                                 int
                                 current_level_modC,
                                 int mC,
                                 const unordered_set<int>&
                                 leaves_levels_modC,
                                 unordered_map<long long, long double>&             cache,
                                 unordered_map<long long, vector <vector<int> > >&  cacheassignment,
                                 bool                                               flexible,
                                 long double                                        GCweight,
                                 int                                                n)
{
  /*
   *  Input:
   * v, a secondary structure tree with a boolean stored in each node for loop assignment 0 <-> proper, 1 <-> unproper
   * c, the assignment for node v[0], 8 possibilities attributes to 8 different int: GC <-> 0, CG <-> 1, AU <-> 2, UA <-> 3, ?? <-> 4, A <-> 5, C <-> 6, RR <-> 7
   * current_level_modA (resp C), the current level modulus mA (resp. mC)
   * mA and mC, the modulo considered for (mA, mC)-separability
   * leaves_levels_modA (resp. C), the set of levels specific to the leaves assigned A (resp. C)
   * cache, the partial solutions for subtrees
   * cacheassignment, the list of local assignments for loops already encountered
   * flexible, a boolean put to 1 if unproper should be the default behavior for every loop
   * GCweight, the weight to attribute to GC base pairs, the default choice should be 0.0
   * n, size of the sequence
   *  Output:
   * A design uniformly sampled for the current tree under these parameters
   */
  long long state = v->parent.first * (n + 1) * 8 * mA * mC + v->parent.second * 8 * mA * mC + c *
                    mA * mC + current_level_modA * mC + current_level_modC;

  if (!cache.count(state)) {
    throw runtime_error("cache miscomputed");
  }

  if (c == 4) {
    long double r = distrib(gen) * cache[state];
    if (v->is_leaf) {
      for (size_t leaf = 5; leaf < 7; leaf++) {
        r -= num_design_biseparated(v,
                                    leaf,
                                    current_level_modA,
                                    mA,
                                    leaves_levels_modA,
                                    current_level_modC,
                                    mC,
                                    leaves_levels_modC,
                                    cache,
                                    cacheassignment,
                                    flexible,
                                    GCweight,
                                    n);
        if (r < 0.0) {
          return stochastic_backtrack_biseparated(v,
                                                  leaf,
                                                  current_level_modA,
                                                  mA,
                                                  leaves_levels_modA,
                                                  current_level_modC,
                                                  mC,
                                                  leaves_levels_modC,
                                                  cache,
                                                  cacheassignment,
                                                  flexible,
                                                  GCweight,
                                                  n);
        }
      }
    } else {
      for (size_t couple = 0; couple < 4; couple++) {
        r -= num_design_biseparated(v,
                                    couple,
                                    current_level_modA,
                                    mA,
                                    leaves_levels_modA,
                                    current_level_modC,
                                    mC,
                                    leaves_levels_modC,
                                    cache,
                                    cacheassignment,
                                    flexible,
                                    GCweight,
                                    n);
        if (r < 0.0) {
          return stochastic_backtrack_biseparated(v,
                                                  couple,
                                                  current_level_modA,
                                                  mA,
                                                  leaves_levels_modA,
                                                  current_level_modC,
                                                  mC,
                                                  leaves_levels_modC,
                                                  cache,
                                                  cacheassignment,
                                                  flexible,
                                                  GCweight,
                                                  n);
        }
      }
    }
  }

  if (v->is_leaf) {
    if (c == 5) {
      return "A";
    }

    if (c == 6) {
      return "C";
    }
  } else {
    int next_level_modA = (mA + (current_level_modA + delta[c]) % mA) % mA;
    int next_level_modC = (mC + (current_level_modC + deltaC[c]) % mC) % mC;
    if ((!(v->unproper)) && (!flexible)) {
      long double r = distrib(gen) * cache[state];
      long long   k = (v->leafbits) * max_degree_optimized_multiloop * 8 + (v->nb_children) * 8 + c;
      if (!cacheassignment.count(k)) {
        throw runtime_error("cache miscomputed");
      }

      const auto& assignments = cacheassignment[k];
      //vector <vector<int>> assignments = cacheassignment[k];
      for (size_t assignment_index = 0; assignment_index < assignments.size(); assignment_index++) {
        long double prod      = 1.0;
        long double GCfactor  = exp(GCweight);
        if ((GCweight != 0.0) && (c == 0 || c == 1)) {
          prod = GCfactor;
        }

        for (size_t i = 0; i < (v->children).size(); i++)
          prod *= num_design_biseparated((v->children)[i],
                                         assignments[assignment_index][i],
                                         next_level_modA,
                                         mA,
                                         leaves_levels_modA,
                                         next_level_modC,
                                         mC,
                                         leaves_levels_modC,
                                         cache,
                                         cacheassignment,
                                         flexible,
                                         GCweight,
                                         n);
        r -= prod;
        if (r < 0.0) {
          string res = "";
          for (size_t i = 0; i < (v->children).size(); i++) {
            string loc = stochastic_backtrack_biseparated((v->children)[i],
                                                          assignments[assignment_index][i],
                                                          next_level_modA,
                                                          mA,
                                                          leaves_levels_modA,
                                                          next_level_modC,
                                                          mC,
                                                          leaves_levels_modC,
                                                          cache,
                                                          cacheassignment,
                                                          flexible,
                                                          GCweight,
                                                          n);
            res += loc;
          }
          if (c == 7) {
            return res;
          } else {
            pair <string, string> p = unwrap(c);
            return p.first + res + p.second;
          }
        }
      }
    } else {
      string res = "";
      for (size_t i = 0; i < (v->children).size(); i++) {
        string loc = stochastic_backtrack_biseparated((v->children)[i],
                                                      4,
                                                      next_level_modA,
                                                      mA,
                                                      leaves_levels_modA,
                                                      next_level_modC,
                                                      mC,
                                                      leaves_levels_modC,
                                                      cache,
                                                      cacheassignment,
                                                      flexible,
                                                      GCweight,
                                                      n);
        res += loc;
      }
      if (c == 7) {
        return res;
      } else {
        pair <string, string> p = unwrap(c);
        return p.first + res + p.second;
      }
    }
  }

  throw invalid_argument("Should not reach here, stochastic backtrack failed");
}


static vector<int>
ssparse(const string& seq)
/*
 *  Input:
 * seq, a secondary structure as a well-parenthesized string
 *  Output:
 * A vector with -1 for an unpaired position and the index of its correspondent for a paired position
 */
{
  vector<int> res(seq.size(), -1);
  stack<int>  p;
  int         j;

  for (size_t i = 0; i < seq.size(); i++) {
    if (seq[i] == '(') {
      p.push(i);
    } else if (seq[i] == ')') {
      j = p.top();
      p.pop();
      res[i]  = j;
      res[j]  = i;
    }
  }
  return res;
}


static bool
checker(const vector<Tree *>& v,
        int                   lim1,
        int                   lim2)
{
  /*
   *  Input:
   * v, a list of secondary structure tree
   * lim1, limits of base pairs for the loops containing no unpaired positions, it should be 3 (and 4 for the root)
   * lim2, limits of base pairs for the loops containing some unpaired positions, it should be 1 (and 2 for the root)
   *  Output:
   * A boolean that is one if the loop is m3o or m5
   */
  int len_BP_children = 0;
  int len_BP_leaves   = 0;

  for (size_t i = 0; i < v.size(); i++) {
    if ((v[i])->is_leaf) {
      len_BP_leaves += 1;
    } else {
      len_BP_children += 1;
    }
  }
  return (len_BP_children > lim1) || (len_BP_leaves > 0 && len_BP_children > lim2);
}


static vector<Tree *>
dbn_to_tree_list(const vector<int>& dbn,
                 int                i,
                 int                j)
/*
 *  Input:
 * dbn, a list with -1 for an unpaired position and the index of its correspondent for a paired position
 * i, j, ints that correspond to the encapsulating base pair that starts the tree (if exists)
 *  Output:
 * A secondary structure tree
 */
{
  vector <Tree *> res;

  while (i <= j) {
    if (dbn[i] == -1) {
      Tree *loc = new Tree;
      loc->parent       = { i, i };
      loc->unproper     = 0;
      loc->children     = vector<Tree *>{};
      loc->leafbits     = 0;
      loc->nb_children  = 0;
      loc->is_leaf      = 1;
      res.emplace_back(loc);
      i += 1;
    } else {
      int             ii    = i;
      int             jj    = dbn[i];
      vector <Tree *> v     = dbn_to_tree_list(dbn, ii + 1, jj - 1);
      Tree            *resu = new Tree;
      resu->parent      = { ii, jj };
      resu->unproper    = checker(v, 3, 1);
      resu->children    = v;
      resu->nb_children = v.size();
      if ((resu->nb_children) > max_degree_optimized_multiloop) {
        resu->leafbits = -1;
      } else {
        resu->leafbits = make_leafbits(v);
      }

      resu->is_leaf = 0;
      res.emplace_back(resu);
      i = dbn[i] + 1;
    }
  }
  return res;
}


static Tree *
dbn_to_tree(const vector<int>& dbn)
/*
 *  Input:
 * dbn, a list with -1 for an unpaired position and the index of its correspondent for a paired position
 *  Output:
 * A secondary structure tree
 */
{
  vector <Tree *> v     = dbn_to_tree_list(dbn, 0, dbn.size() - 1);
  Tree            *resu = new Tree;

  resu->parent      = { 0, dbn.size() };
  resu->unproper    = checker(v, 4, 2);
  resu->children    = v;
  resu->nb_children = v.size();
  if ((resu->nb_children) > max_degree_optimized_multiloop) {
    resu->leafbits = -1;
  } else {
    resu->leafbits = make_leafbits(v);
  }

  resu->is_leaf = 0;
  return resu;
}


static vector<vector<int> >
part(int i)
/*
 *  Input:
 * i, an int
 *  Output:
 * Enumerate all sets of i + 1 elements
 */
{
  if (i == -1) {
    return { {} };
  }

  vector<vector<int> >  li = part(i - 1);

  vector<vector<int> >  res = li;

  for (vector<int> elem : li) {
    elem.push_back(i);
    res.push_back(elem);
  }
  return res;
}


static void
tree_cleanup(Tree *v)
{
  if (!v) {
    return;
  }

  for (Tree *t:v->children)
    tree_cleanup(t);
  delete v;
}


static unordered_set<string>
sample(size_t                                             magic_threshold,
       long double                                        Zfull,
       const Tree                                         *tree,
       int                                                m,
       vector<unordered_set<int> >                        tuple_leaf_level_list,
       vector<long double>                                Zlist,
       vector<unordered_map<long long, long double> >     cachelist,
       unordered_map<long long, vector <vector<int> > >&  cacheassignment,
       bool                                               flex,
       long double                                        GCweight,
       int                                                n)
/*
 *  Input:
 * magic_threshold, the number of time to sample to get sequences
 * Zfull, the sum of partition functions at a given modulus m
 * tree, a secondary structure tree with a boolean stored in each node for loop assignment 0 <-> proper, 1 <-> unproper
 * tuple_leaf_level_list, the list of set of leaf levels
 * Zlist, the list of partition function for each set of leaf levels
 * m, the modulus considered for m-separability
 * cachelist, the partial solutions for subtrees for each set of leaf levels
 * cacheassignment, the list of local assignments for loops already encountered
 * flex, a boolean put to 1 if unproper should be the default behavior for every loop
 * GCweight, the weight to attribute to GC base pairs, the default choice should be 0.0
 * n, size of the sequence
 *  Output:
 * A set of designs of size at most magic_threshold uniformly sampled at level m for the current tree under these parameters
 */
{
  unordered_set<string> resu;

  if (Zfull == 0.0) {
    return resu;
  }

  for (size_t it = 0; it < magic_threshold; it++) {
    long double r = distrib(gen) * Zfull;
    for (size_t i = 0; i < Zlist.size(); i++) {
      r -= Zlist[i];
      if (r < 0.0) {
        resu.insert(stochastic_backtrack(tree,
                                         7,
                                         0,
                                         m,
                                         tuple_leaf_level_list[i],
                                         cachelist[i],
                                         cacheassignment,
                                         flex,
                                         GCweight,
                                         n));
        break;
      }
    }
  }
  return resu;
}


static unordered_set<string>
sample_biseparated(size_t                                             magic_threshold,
                   long double                                        Zfull,
                   const Tree                                         *tree,
                   int                                                mA,
                   int                                                mC,
                   vector<unordered_set<int> >                        tuple_leaf_levelA_list,
                   vector<unordered_set<int> >                        tuple_leaf_levelC_list,
                   vector<long double>                                Zlist,
                   vector<unordered_map<long long, long double> >     cachelist,
                   unordered_map<long long, vector <vector<int> > >&  cacheassignment,
                   bool                                               flex,
                   long double                                        GCweight,
                   int                                                n)
/*
 *  Input:
 * magic_threshold, the number of time to sample to get sequences
 * Zfull, the sum of partition functions at given modulo (mA, mC)
 * tree, a secondary structure tree with a boolean stored in each node for loop assignment 0 <-> proper, 1 <-> unproper
 * tuple_leaf_levelA_list/tuple_leaf_levelC_list, the list of set of leaf levels for A-levels/C-levels
 * Zlist, the list of partition function for each set of leaf levels
 * (mA, mC), the modulo considered for (mA, mC)-separability
 * cachelist, the partial solutions for subtrees for each set of leaf levels
 * cacheassignment, the list of local assignments for loops already encountered
 * flex, a boolean put to 1 if unproper should be the default behavior for every loop
 * GCweight, the weight to attribute to GC base pairs, the default choice should be 0.0
 * n, size of the sequence
 *  Output:
 * A set of designs of size at most magic_threshold uniformly sampled at level  (mA, mC) for the current tree under these parameters
 */
{
  unordered_set<string> resu;

  if (Zfull == 0.0) {
    return resu;
  }

  for (size_t it = 0; it < magic_threshold; it++) {
    long double r = distrib(gen) * Zfull;
    for (size_t i = 0; i < Zlist.size(); i++) {
      r -= Zlist[i];
      if (r < 0.0) {
        resu.insert(stochastic_backtrack_biseparated(tree,
                                                     7,
                                                     0,
                                                     mA,
                                                     tuple_leaf_levelA_list[i],
                                                     0,
                                                     mC,
                                                     tuple_leaf_levelC_list[i],
                                                     cachelist[i],
                                                     cacheassignment,
                                                     flex,
                                                     GCweight,
                                                     n));
        break;
      }
    }
  }
  return resu;
}


vector<string>
generate_seeds(string                       structure,
               unsigned int                 target_nb_sols,
               vrna_inverse_linearbp_opt_t  *p)
{
  if (target_nb_sols == 0) {
    vector<string> resu;
    return resu;
  }

  int         n               = structure.size();
  size_t      magic_threshold = 10 * target_nb_sols * log(target_nb_sols);  //Coupon collector times 10 for safety
  bool        biseparable     = 0;
  int         maxmodulo       = 3;                                          //Maximal modulos considered for separability or biseparability
  long double GCweight        = 0.0;
  if (p) {
    biseparable = p->biseparable;
    maxmodulo   = p->maxmodulo;
    GCweight    = p->GCweight;
  }

  Tree                                              *tree = dbn_to_tree(ssparse(structure));
  unordered_map<long long, long double>             cache;
  unordered_map<long long, vector <vector<int> > >  cacheassignment;
  unordered_set <string>                            preresu;
  bool                                              done = 0;
  for (int flex = 0; flex < 2; flex++) {
    cacheassignment.clear();
    if (biseparable) {
      int                 mA;
      int                 mC;
      unordered_set<int>  tuple_leaf_levelA;
      unordered_set<int>  tuple_leaf_levelC;

      for (int sumM = 4; sumM < 2 * maxmodulo + 1; sumM++) {
        for (mA = 2; mA < sumM - 1; mA++) {
          mC = sumM - mA;
          vector<long double>                             Zlist;
          vector<unordered_map<long long, long double> >  cachelist;
          vector<unordered_set<int> >                     tuple_leaf_levelA_list;
          vector<unordered_set<int> >                     tuple_leaf_levelC_list;
          vector<vector<int> >                            liA   = part(mA - 1);
          vector<vector<int> >                            liC   = part(mC - 1);
          long double                                     Zfull = 0.0;
          for (vector<int> leaf_levelA : liA) {
            for (vector<int> leaf_levelC : liC) {
              cache.clear();
              tuple_leaf_levelA = unordered_set<int>(leaf_levelA.begin(), leaf_levelA.end());
              tuple_leaf_levelC = unordered_set<int>(leaf_levelC.begin(), leaf_levelC.end());
              long double Z = num_design_biseparated(tree,
                                                     7,
                                                     0,
                                                     mA,
                                                     tuple_leaf_levelA,
                                                     0,
                                                     mC,
                                                     tuple_leaf_levelC,
                                                     cache,
                                                     cacheassignment,
                                                     flex,
                                                     GCweight,
                                                     n);
              Zfull += Z;
              if (Z > 0.0) {
                Zlist.push_back(Z);
                cachelist.push_back(cache);
                tuple_leaf_levelA_list.emplace_back(tuple_leaf_levelA);
                tuple_leaf_levelC_list.emplace_back(tuple_leaf_levelC);
              }
            }
          }
          preresu.merge(sample_biseparated(magic_threshold,
                                           Zfull,
                                           tree,
                                           mA,
                                           mC,
                                           tuple_leaf_levelA_list,
                                           tuple_leaf_levelC_list,
                                           Zlist,
                                           cachelist,
                                           cacheassignment,
                                           flex,
                                           GCweight,
                                           n));
          if (preresu.size() >= target_nb_sols) {
            done = 1;
          }

          if (done) {
            break;
          }
        }
        if (done) {
          break;
        }
      }
    } else {
      int                 m;
      unordered_set<int>  tuple_leaf_level;
      for (m = 2; m < maxmodulo + 1; m++) {
        vector<long double>                             Zlist;
        vector<unordered_map<long long, long double> >  cachelist;
        vector<unordered_set<int> >                     tuple_leaf_level_list;
        vector<vector<int> >                            li    = part(m - 1);
        long double                                     Zfull = 0.0;
        for (vector<int> leaf_level : li) {
          cache.clear();
          tuple_leaf_level = unordered_set<int>(leaf_level.begin(), leaf_level.end());
          long double Z = num_design(tree,
                                     7,
                                     0,
                                     m,
                                     tuple_leaf_level,
                                     cache,
                                     cacheassignment,
                                     flex,
                                     GCweight,
                                     n);
          Zfull += Z;
          if (Z > 0.0) {
            Zlist.push_back(Z);
            cachelist.push_back(cache);
            tuple_leaf_level_list.emplace_back(tuple_leaf_level);
          }
        }
        preresu.merge(sample(magic_threshold,
                             Zfull,
                             tree,
                             m,
                             tuple_leaf_level_list,
                             Zlist,
                             cachelist,
                             cacheassignment,
                             flex,
                             GCweight,
                             n));
        if (preresu.size() >= target_nb_sols) {
          done = 1;
        }

        if (done) {
          break;
        }
      }
    }
  }
  vector<string> resu;
  resu.insert(resu.end(), preresu.begin(), preresu.end());
  resu.resize(target_nb_sols);
  tree_cleanup(tree);
  return resu;
}
