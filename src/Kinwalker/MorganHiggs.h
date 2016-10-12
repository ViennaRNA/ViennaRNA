/*
  Last changed Time-stamp: <2006-02-14 01:01:41 xtof>
  $Id: MorganHiggs.h,v 1.7 2007/10/01 08:29:47 Kinwalker Exp $
*/
#ifndef _MORGANHIGGS_H_
#define _MORGANHIGGS_H_
#include <deque>
#include <vector>
#include <string>
#include <utility>
#include <map>
#include <algorithm>
#include <iterator>

#include <cmath>
#include <limits>
#include <iostream>
#include <fstream>

#include "Energy.h"
#include "Util.h"
#include "options.h"

//std::pair<double,std::string>
std::vector<std::pair<double,std::string> >
MorganHiggsStudlaEnergy(std::string sequence,
		  std::string src,
		  std::string tgt,
		  double saddlE,int look_ahead,std::string grouping);

std::vector<std::pair<double,std::string> >
MorganHiggsEnergy(std::string sequence,
		  std::string src,
		  std::string tgt,
		  double saddlE,int look_ahead,std::string grouping);

//std::map<int,std::vector<int> > GetConflictGroup(std::vector<std::pair<int,int> > only_in_base_pairs,std::vector<std::pair<int,int> > only_in_base_pairs2);
void GetConflictGroup(std::map<int,std::vector<int> > & conflict_group, const std::vector<std::pair<int,int> > & only_in_base_pairs,const std::vector<std::pair<int,int> > & only_in_base_pairs2);

//std::vector<std::pair<double,std::string> > PartialPath(std::vector<int> combination, std::string sequence,std::vector<std::pair<int,int> > backtrack_base, 
//						std::map<int,std::vector<int> > conflict_group,std::vector<std::pair<int,int> > only_in_base_pairs,std::vector<std::pair<int,int> > only_in_base_pairs2);

void DoPartialPath(std::vector<std::pair<double,std::string> > & path, const std::vector<int> & combination, std::string sequence, const std::vector<std::pair<int,int> > & 
backtrack_base,  const std::map<int,std::vector<int> > & conflict_group, const std::vector<std::pair<int,int> > & only_in_base_pairs,
		  const std::vector<std::pair<int,int> > & only_in_base_pairs2);

//void PartialPath(std::vector<std::pair<double,std::string> > & path,const std::vector<int> & combination, std::string sequence,const std::vector<std::pair<int,int> > & 
//backtrack_base, const std::map<int,std::vector<int> > & conflict_group,const std::vector<std::pair<int,int> > & only_in_base_pairs,
//		 const std::vector<std::pair<int,int> > & only_in_base_pairs2);


//std::string PrintConflictGroup(std::map<int,std::vector<int> > conflict_group);
std::string PrintConflictGroup(std::map<int,std::vector<int> > conflict_group,std::vector<std::pair<int,int> > only_in_base_pairs,
			       std::vector<std::pair<int,int> > only_in_base_pairs2);

std::string PrintConflictGroup(std::vector<std::vector<int> > conflict_group,std::vector<std::pair<int,int> > only_in_base_pairs,
			       std::vector<std::pair<int,int> > only_in_base_pairs2);

std::string PrintCombination(std::vector<int> combination);
void MakePairTableFromBasePairs(const std::vector<std::pair<int,int> >& base_pairs,int length);

 void UpdateBacktrackBase(std::vector<std::pair<int,int> >& base_pairs,const std::map<int,std::vector<int> > & conflict_group,const std::vector<int> & remove,
			  const std::vector<std::pair<int,int> > & only_in_base_pairs,const std::vector<std::pair<int,int> > & only_in_base_pairs2);

void FindBestPartialPathCombination(std::vector<int> & best_combination,int lookahead,std::string sequence,const std::vector<std::pair<int,int> > & backtrack_base,
const std::map<int,std::vector<int> > & conflict_group,const std::vector<std::pair<int,int> > & only_in_base_pairs,
				    const std::vector<std::pair<int,int> > & only_in_base_pairs2);

bool IsGCPair(std::pair<int,int> bp,const std::string & sequence);

std::vector<std::pair<int,int> > GetLonelyBasePairs(const std::vector<std::pair<int,int> > & current,const std::vector<std::pair<int,int> > & target,const std::string & sequence );









/**
class for a base pair stack, i.e. sets of consecutive basepairs. Bulges do not
count.
*/
class Stack{
public:

  
  //are all bases in the stack GC
  bool GC;//=true;

  //maximal size for which a GC stack is minuscule
  const static int GCMaxMinusculeStackSize=1;


  //maximal size for which a non-GC stack is minuscule
  const static int NonGCMaxMinusculeStackSize=2;
 

  //yhr base pairs that make up the stack. the innermost bp is at the front, the outermost at the back.
  std::deque<std::pair<int,int> > bp;

  
  bool CanBeAddedInside(std::pair<int,int> bp );

  bool CanBeAddedOutside(std::pair<int,int> bp );

  bool CanBeAddedInside(int first,int second);

  bool CanBeAddedOutside(int first,int second);

  bool IsMinuscule();

  static bool IsMinuscule(int size,bool GC);

  Stack();

 
  int GetSize();


 
  bool Contains(std::pair<int,int> bp );

  bool Contains(int first,int second);

  void Remove(std::pair<int,int> bp );

  void Remove(int bp_idx);


  void AddOutside(std::pair<int,int> bp,char first,char second);

  void AddInside(std::pair<int,int> bp,char first,char second);

  void AddOutside(int first,int second,char cfirst,char csecond);

  void AddInside(int first,int second,char cfirst,char csecond);


 static  bool IsGCPair(char first,char second);

  void Clear();

  

  static inline bool IsGC(char base);

  std::string  Print();
};








/**
PartialPath tries all combination in a conflict group and appends the best path to the path traversed so far
NOTE:Indices of elements that are added are augmented by BP_ADD_CONST, so when undoing a trail, we now in which
data structure the index was and whether it is to be added or removed.
*/

class PartialPath{
public:
 
 PartialPath(std::vector<std::pair<int,int> > src,std::vector<std::pair<int,int> > target,std::string sequence,int lookahead,const std::vector<std::pair<int,int> > & only_in_base_pairs,const std::vector<std::pair<int,int> > & only_in_base_pairs2);


 
  /** 
    Function invoked to find partial paths and append them to sequence of structures. Creates a conflict group, then tries all combinations and selects a 
     path covering lookahead members of the conflict group until it is empty.
  */

  void FindOptimalPath(std::vector<std::pair<double,std::string> >& structures);

  /**
      Add/removes all elements in add/remove_catalog given in combination.
  */

  double WalkPartialPath( std::vector<int> combination );

  /**
     adds/removes all elements in the conflit group at the given index thus constructing
a partial path
  */
  void TreatConflictGroupIndex(int idx);


  /*
   returns true if the pair (first,second) is in pair_table, false otherwise.
  */

  bool PairIsInPairtable(int first,int second);


  /**
     Get the first element of the base pair at index idx in the add_catalog,
 i.e. those elts in the conflict group that need ot be added to the src structure
  */
  // int GetFirstFromAddCatalog(int idx);
  /**
     Get the second element of the base pair at index idx in the add_catalog,
  */
  //  int GetSecondFromAddCatalog(int idx);

  /*
     Retrieves index of the first element in the idx-th component of conflict_group,
     which is an index in the add_catalog
  */
 int GetAddCatalogIndex(int idx);


 /*
     Retrieves index of the remove_idx element in the group_idx-th component of conflict_group,
     which is an index in the remove_catalog
  */
  int GetRemoveCatalogIndex(int group_idx,int remove_idx);
 

 /**
     Get the first element at index of the remove_idx element in the group_idx-th component of conflict_group
     i.e. remove_catalog[conflict_group[group_idx][remove_idx]].first
  */

  int GetFirstFromConflictGroupRemove(int group_idx,int remove_idx);


 /**
     Get the second element at index of the remove_idx element in the group_idx-th component of conflict_group
     i.e. remove_catalog[conflict_group[group_idx][remove_idx]].second
  */

  int GetSecondFromConflictGroupRemove(int group_idx,int remove_idx);

  /*
     Retrieves first element of base pair in add_catalog at index of the first element in the idx-th component of conflict_group,
     i.e. add_catalog[conflict_group[group_idx][0]].first
  */

 int GetFirstFromConflictGroupAdd(int group_idx);


  /*
     Retrieves second element of base pair in add_catalog at index of the first element in the idx-th component of conflict_group,
     i.e. add_catalog[conflict_group[group_idx][0]].first
  */

  int GetSecondFromConflictGroupAdd(int group_idx);


 /**
     add the base pair bp to pair_table
  */
  //  void AddBasePair(const std::pair<int,int> bp);

  /**
    Add a pair to pair_table
  */
  void AddToPairTable(std::pair<int,int> bp);

  /**
    Add the pair (first,second) to pair_table.
  */
  void AddToPairTable(int first,int second);

 /**
    Add the pair (first,second) to pair_table and record to this->currentPartialPathTrail so it can be undone by UndoPartialPathTrail().
  */
  void AddToPairTableWithLog(int first,int second,int addCatalogIdx);


  /**
    Add pair from pair_table
  */

  void RemoveFromPairTable(std::pair<int,int> bp);
  /**
    Remove the pair (first,second) from pair_table.
  */
  void RemoveFromPairTable(int first,int second);

/**
    Remove the pair (first,second) to pair_table and record to this->currentPartialPathTrail so it can be undone by UndoPartialPathTrail().
  */
  void RemoveFromPairTableWithLog(int first,int second,int removeCatalogIdx);

  /**
    Lists all minuscule Stacks at the beginning of treating a conflict group. This ensures that after the removal
    of initial minuscule stacks there is at most 1 minuscule stack at a time.
  */
  std::vector<Stack> GetMinusculeStacks();

  //=========================== 
  // list of base_pairs that need to be added to move from src to target. They make up the first index in each component of a conflict_group
  std::vector<std::pair<int,int> > add_catalog;
 // list of base_pairs that need to be removed to move from src to target. They make up all indices after the first in each component of a conflict_group
  std::vector<std::pair<int,int> > remove_catalog;


  /*
   List of vectors. The first element in each vector gives index of element in add_catalog, all other indices point to elts in remove_catalog that are
   in conflict with it.
  */
  std::vector<std::vector<int> > conflict_group;


  /**
    A sequence of this->lookahead indices to elements in conflict_group that denote a partial path.
  */
  std::vector<int> combination;

  
  /**
    The initial structure for the barrier heuristic
  */
  std::vector<std::pair<int,int> > src;
 
  /*
  final structure of barrier heuristic
  */
 std::vector<std::pair<int,int> > target;

  /*
     The RNA that is being folded
  */

     std::string sequence;

  /**
    Number of elements treated at once in a combination.
  */
  int lookahead;
 

  /**
    the current minuscule stack (at most one is allowed). If there
is none, it is empty, i.e. of size 0.
  */
 
  Stack minusculeStack;

 
  /**
     list of indices add/removed during latest call to TreatConflictGroupIndex()
  */
  std::vector<int> currentPartialPathTrail;

  /**
     list of indices add/removed during the call to TreatConflictGroupIndex() that yielded the
     lowest saddle of all combinations tried so far.
  */

  std::vector<int> bestPartialPathTrail;

  /**
     list of concatenated partial paths since the conflict_group has been built. (i.e. FindOptimalPath was called the last time). 
  */
  std::vector<int> pathTrail;

 


  /**
     Add the structures traversed by integrating the elements from conflict_group
     to the list of such structures. called at the end of FindOptimalPath.
  */

  void AddPathTrailToStructureList(std::vector<std::pair<double,std::string> >& structures);
 
  /**
     Add the indices of the elts in add/remove_catalog that were added to the path when dealing
     with the best combination of lookahead elements in conflict_group to the path trail.
     Updates pair_table to the structure obtained when best partial path is taken.
  */
  void  AddBestPartialPathTrailToPathTrail();

  /**
     Empty the vector that holds the last combination of treated indices for a combination
     of lookahead elelemtns from the conflict_group.
  */
  void UndoPartialPathTrail();

  
  /*
     Find the elements in add_catalog that have the lowest number of conflicts with the current
     structure and build the conflict_group form those and their conflicts.

  */
  void GetConflictGroup();

  /*
    Remove the base pairs in the minuscule Stack mStack from pair_table unless it is part of the target conformation.
    In that case base_pairs are added (inside first) until it is not minuscule any more.

  */

  void RemoveMinusculeStack(Stack mStack);
  /*
    Remove the base pairs in the instnace variable mStack from pair_table unless it is part of the target conformation.
    In that case base_pairs are added (inside first) until it is not minuscule any more.
    
  */

  void RemoveMinusculeStack();


  /**
     Returns true if a adding the base pair (first,second) creates a second minuscule stack, false otherwise.
     Note that adding a base pair can never create two minuscle stacks. 
  */
  bool BPAdditionCreatesSecondMinusculeStack(int first,int second);

  /**
    Update add_catalog, remove_catalog and src after the best path for the elements in conflict_group has been
    found and accepted as a path.
   */
  void UpdateBacktrackBase();

  /**
    Counts the number of stacks that are created by removing the base pair (first,second) from pair_table.
    If a new minuscule stakc is created it is placed in createdMinusculeStack1, a second one in createdMinusculeStack2
    (both global variables). If mStack disappears by removing (first,second) it is cleared.
    
  */

  int StacksCreatedByBPRemoval(int first,int second);

  /**
  Returns true if the pair (first,second) is lonely in pair_table, false otherwise.
  */
  bool IsLonelyPair(int first,int second);

  /**
   Returns true if the minuscul stack is contained in target. Actually, if both the
  first and last bp are in target, the function returns true as currently no stacks
  greater than size 2 are minuscule. Returns false if the stack is not contained in the
  target structure.
  */
  bool MinusculeStackIsInTarget(Stack mStack);
};







#endif

/* End of file */
