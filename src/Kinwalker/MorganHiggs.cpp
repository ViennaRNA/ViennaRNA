/*
  Last changed Time-stamp: <2007-04-02 14:38:35 caro>
  $Id: MorganHiggs.cpp,v 1.13 2007/10/21 21:01:35 Kinwalker Exp $
*/

#include <cstring>
#include "MorganHiggs.h"
#define BP_ADD_CONST 10000
using std::cout;
using std::endl;

extern short * pair_table;

extern "C" {
#include "energy_const.h" /* defines INF */
}

bool MHS_debug=false;
int min_stack_size=1;

/**
  rekursive Fakultaetsformel fuer n>=0.
*/
inline int
Factorial(int n)
{
  if(n==1 || n==0)
    return (1);

  return (n*Factorial(n-1));
}

int N_take_k(int n,int k){

  if(k>=n) return Factorial(n);
  int ret=1;
  for(int i=0;i<k;i++) ret*=n-i;
  return ret;

}


/**
 * Indices in the permutation group start at 0
*/
std::vector<int>
GetPermutation(int elts,int idx)
{
  std::vector<int> permutation = std::vector<int>();

  for (int i=0; i<elts-1; i++) {
    
    // position i is occupied by element rank, i.e. which takes the
    // slice of element [(elts-i)*rank,(elts-i)*(rank+1)] in the
    // lexicographical order of the elts in the permutation group
    int basket_size = Factorial(elts-i-1);
    int rank = (int)std::floor(idx/basket_size);
    permutation.push_back(rank);
    idx = idx % basket_size;
  }
  
  std::vector<int> ret = std::vector<int>();
  std::vector<int> available = std::vector<int>();

  for (int i=1; i<=elts; i++) available.push_back(i);

  for(int i=0;i<elts-1;i++) {
    // have to choose the permutation[i]'s elements from those that
    // are still available.
    ret.push_back(available[permutation[i]]);
    available.erase(available.begin()+permutation[i]);
  }
  
  // with the last element there is no choice.
  ret.push_back(available.front());
  return ret;
}


/**
  Get a combination of k elements from a set of n without replacement.
  Idx denotes which of the n choose k (=binomial coefficient) combinations is
  returned. Indices range from 1 to n.
*/

std::vector<int> GetCombination(int n,int k,long int idx)
{
  
  if(n<=k) return GetPermutation(n,idx);
  #ifdef _DEBUG_MH_GETCOMBINATION_
      std::cout<<"n "<<n<<"\n";
      std::cout<<"k "<<k<<"\n";
      std::cout<<"idx "<<idx<<"\n";
  #endif
  //Example 10 choose 4
  //idx=x_0*7*8*9+x_1*7*8+x_2*7+x_3
  //9*8*7
  long int modulus=1;
  for(int i=n-1;i>n-k;i--) modulus*=i;

  #ifdef _DEBUG_MH_GETCOMBINATION_  
      std::cout<<"final modulus "<<modulus<<"\n";
  #endif
  // vector<int> ret(n);

  //std::vector<int> elts=std::vector<int>(); 
  //workaround as the the previous line did not compile
  //std::vector<int> elts = *(new std::vector<int>());
      std::vector<int> elts = std::vector<int>();
  for (int i=1; i<=n; i++) elts.push_back(i);
  //i is idx of a number in the combination, i=0 stands for the right side.
  for(int i=0;i<k;i++){
    #ifdef _DEBUG_MH_GETCOMBINATION_  
        std::cout<<"i "<<i<<"\n";
    #endif
    //ranges from 0 to n-i-1
    #ifdef _DEBUG_MH_GETCOMBINATION_  
	std::cout<<"idx "<<idx<<"\n";
    #endif
    int chosen_idx=(int)floor(idx/modulus);
    #ifdef _DEBUG_MH_GETCOMBINATION_  
        std::cout<<"chosen_idx "<<chosen_idx<<"\n";
    #endif
    idx-=chosen_idx*modulus;
    #ifdef _DEBUG_MH_GETCOMBINATION_  
        std::cout<<"idx "<<idx<<"\n";
    #endif
    //if condition to avoid division by 0
    if(modulus!=1) modulus=(int)(modulus/(n-i-1));
    #ifdef _DEBUG_MH_GETCOMBINATION_ 
        std::cout<<"modulus "<<modulus<<"\n";
    #endif
    int i1=elts[i];
    int i2=elts[i+chosen_idx];
    #ifdef _DEBUG_MH_GETCOMBINATION_ 
        std::cout<<"i1 "<<i1<<"\n";
    #endif
    #ifdef _DEBUG_MH_GETCOMBINATION_ 
	std::cout<<"i2 "<<i2<<"\n";
    #endif
    elts[i]=i2;
    elts[i+chosen_idx]=i1;
  }
  elts.erase(elts.begin()+k,elts.end());
  #ifdef _DEBUG_MH_GETCOMBINATION_ 
  std::cout<<"GetCombination returns"<<std::endl;
  #endif
  return elts;
}
/**
  Returns true, if p1 and p2 are incompatible, false otherwise.
 Returns false if p1==p2.
*/

bool BasePairConflict(std::pair<int,int> p1, std::pair<int,int> p2){
   // same base pair in both structures
    // (-----)
    // (-----)
    if (p1.first == p2.first && p1.second == p2.second)
      return (false);

    // base pairs have on position in common
    // (-----)                (-----)
    //       (-----) or (-----)      
    if (p1.first == p2.first || p1.first == p2.second
	|| p1.second == p2.first || p1.second == p2.second)
      return (true);

    // "crossing" base pairs
    // (-------)
    //     (-----)
    if (p1.first <= p2.first && p1.second >= p2.first
	&& p1.second <= p2.second)
      return (true);

    // "crossing" base pairs
    //     (-----)    (---------)
    // (-------)   or (-------)
    if (p1.first >= p2.first && p1.first <= p2.second
	&& p1.second >= p2.second)
      return (true);
    else
      return (false);
}


/**
Lists the indices in a vector pairs that are in conflict with a base pair p.
*/

std::vector<int> 
FindConflicts(std::pair<int,int> p, std::vector<std::pair<int,int> > pairs)
{
  std::vector<int> conflicts= std::vector<int>();

  for (size_t i=0; i<pairs.size(); i++) {
    
    if ( BasePairConflict(p,pairs[i]))
      conflicts.push_back(i);
  }
  
  return (conflicts);
}



/*
  Given a base pair and a sequence, returns if the pair is GC or CG,
  false otherwise. Used by the stack class to determine whether a stack
  i minuscule, as the max size for minuscule GC stacks is Stack::GCMaxMinusculeStackSize
  while it is Stack::NonGCMaxMinusculeStackSize for non-GC stacks.
*/

bool IsGCPair(std::pair<int,int> bp,const std::string & sequence){

  bool first=(sequence[bp.first-1]=='G' || sequence[bp.first-1]=='C');
  bool second=(sequence[bp.second-1]=='G' || sequence[bp.second-1]=='C');

  if(first && second) return true;
  else return false;
}



//default constructor for a "stack" without base pairs
Stack::Stack(){
   GC=true;
}

//removes all basepairs from the stack, so its size becomes 0. GC therefore true.
void Stack::Clear(){
  bp.clear();
  GC=true;
}

//returns number of base pairs in the stack
int Stack::GetSize(){
  return bp.size();
}


/*
  returns true if the base pair (first,second) is inside of and adjacent to the current innermost pair,
false otherwise.
*/
bool Stack::CanBeAddedInside(int first,int second){
  if(GetSize()==0) return true;
  //have to go back from the 5' end if bp is inside innerMost
  if(first-1==bp[0].first && second+1==bp[0].second) return true;
  else return false;
}

/*
  returns true if the base pair (first,second) is outside of and adjacent to the current outermost pair,
false otherwise.
*/
bool Stack::CanBeAddedOutside(int first,int second){
  if(GetSize()==0) return true;
  if(first+1==bp[0].first && second-1==bp[0].second) return true;
  else return false;
}

/*
  returns true if base_pair is inside of and adjacent to the current innermost pair,
false otherwise.
*/
bool Stack::CanBeAddedInside(std::pair<int,int> base_pair ){
  //have to go back from the 5' end if bp is inside innerMost
  if(std::make_pair(base_pair.first-1,base_pair.second+1)==bp[0]) return true;
  else return false;
}

/*
  returns true if base_pair is outside of and adjacent to the current outermost pair,
false otherwise.
*/
bool Stack::CanBeAddedOutside(std::pair<int,int> base_pair ){
  //have to go forward from the 5' end if bp is outside outerMost
  if(std::make_pair(base_pair.first+1,base_pair.second-1)==bp.back()) return true;
  else return false;
}


/**
  ToString() function that lists whether the stack is GC and its basepairs.
*/
std::string  Stack::Print(){
  std::string s;
  s+="GC "+Str(GC);
  for(size_t i=0;i<bp.size();i++) s+=" "+PrintBasePair(bp[i]);
    return s+"\n";
}

/*
 returns true, if all pairs are GC pairs, or there are no pairs.
 False otherwise.
*/
bool Stack::IsGC(char base){
  if(base=='G' || base=='C') return true;
  else return false;
}
/*
Returns true, if first and second are both either 'G' or 'C',false otherwise.
*/
 
bool Stack::IsGCPair(char first,char second){
  if(IsGC(first) && IsGC(second)) return true;
else return false;

}


/*
Returns true, if the stack is minuscule, i.e. of size at most
GCMaxMinusculeStackSize for GC stacks and NonGCMaxMinusculeStackSize 
for non-GC stacks.Static version.
*/

bool Stack::IsMinuscule(int size,bool GC){
  if(GC){
    if(size<=GCMaxMinusculeStackSize) return true;
  }
  else if(size<=NonGCMaxMinusculeStackSize) return true;
  return false;
}

/*
Returns true, if the stack is minuscule, i.e. of size at most
GCMaxMinusculeStackSize for GC stacks and NonGCMaxMinusculeStackSize 
for non-GC stacks.
*/
bool Stack::IsMinuscule(){
  int size=GetSize();
  if(GC){
    if(size<=GCMaxMinusculeStackSize) return true;
  }
  else if(size<=NonGCMaxMinusculeStackSize) return true;
  return false;

}


/*
Returns true, if the base pair (first,second) is part of the stack, false otherwise.
*/
bool Stack::Contains(int first,int second){
   if(find(bp.begin(),bp.end(),std::make_pair(first,second))!=bp.end()) return true;
   else return false;
}

/*
Returns true, if base_pair is part of the stack, false otherwise.
*/
bool Stack::Contains(std::pair<int,int> base_pair ){
   if(find(bp.begin(),bp.end(),base_pair)!=bp.end()) return true;
   else return false;
}

/**
removes base_pair from the stack (provided it was part of it).
*/
void Stack::Remove(std::pair<int,int> base_pair ){
  std::deque<std::pair<int,int> >::iterator it=find(bp.begin(),bp.end(),base_pair);
  bp.erase(it);
}

/**
removes the  base pair at index bp_index in bp from the stack. No protection for invalid
indices.
*/
void Stack::Remove(int bp_index){
  bp.erase(bp.begin()+bp_index);
}



/**
adds the pair base_pair with characters first and second to the outside of the
stack. You need to be sure that the pair actually can be added before calling this
function.
*/

void Stack::AddOutside(std::pair<int,int> base_pair,char first,char second ){
    bp.push_back(base_pair);
    if(!IsGCPair(first,second)) GC=false;
  }

/**
adds the pair (first,second) with characters cfirst and csecond to the outside of the
stack. You need to be sure that the pair actually can be added before calling this
function.
*/
void Stack::AddOutside(int first,int second,char cfirst,char csecond ){
    bp.push_back(std::make_pair(first,second));
    if(!IsGCPair(cfirst,csecond)) GC=false;
  }

/**
adds the pair base_pair with characters first and second to the inside of the
stack. You need to be sure that the pair actually can be added before calling this
function.
*/
  void Stack::AddInside(std::pair<int,int> base_pair,char first,char second ){
    bp.push_front(base_pair);
    if(!IsGCPair(first,second)) GC=false;
  }

/**
adds the pair (first,second) with characters cfirst and csecond to the inside of the
stack. You need to be sure that the pair actually can be added before calling this
function.
*/
 void Stack::AddInside(int first,int second,char cfirst,char csecond ){
    bp.push_front(std::make_pair(first,second));
    if(!IsGCPair(cfirst,csecond)) GC=false;
  }



  

/**
larger one of the potentially two creatd minuscule stacks form base pair removal
*/

Stack createdMinusculeStack1;

Stack createdMinusculeStack2;



void PartialPath::GetConflictGroup() {
  conflict_group.clear();
  //bool debug=false;//true;
  std::vector<int> conflicts=std::vector<int>(add_catalog.size()); 
  // std::map<int,std::vector<int> >  conflict_group=std::map<int,std::vector<int> >();
  if(add_catalog.size()>0){
    size_t minsize=remove_catalog.size();
    for (size_t i=0; i<add_catalog.size();i++) {
      conflicts[i]=0;
      for (size_t j=0; j<remove_catalog.size();j++) {
        if(BasePairConflict(add_catalog[i],remove_catalog[j])) conflicts[i]++;
      }
      #ifdef _DEBUG_MH_
          std::cout<<"add_catalog["+Str((int)i)+"]: ("+Str(add_catalog[i].first)+","+Str(add_catalog[i].second)+")\n";
          std::cout<<"conflicts[i]:"<<conflicts[i]<<"\n";
      #endif
	  if(conflicts[i]<(int)minsize) minsize=conflicts[i]; 
    }
    #ifdef _DEBUG_MH_
        std::cout<<"minsize:"<<minsize<<"\n";
    #endif
    for (size_t i=0; i<add_catalog.size();i++) {
      if(conflicts[i]==(int)minsize){
        std::vector<int> new_conflicts=std::vector<int>();
        //idx in oibp2 (add_catalog) of the base pair to be added 
	new_conflicts.push_back(i);
        //conflict_group[i]=std::vector<int>(); 


        for (size_t j=0; j<remove_catalog.size();j++) {
          if(BasePairConflict(add_catalog[i],remove_catalog[j])) 
             //idx in oibp (remove_catalog) of the base pair to be removed 
     	     new_conflicts.push_back(j);
	    //conflict_group[i].push_back(j);
        }
	conflict_group.push_back(new_conflicts);
      }
    }
  }
  //all has been added, so only bp to be removed are left
  else{
    //conflict_group[-1]=std::vector<int>();
    std::vector<int> new_conflicts=std::vector<int>();
    //idx in oibp2 (add_catalog) of the base pair to be added 
    new_conflicts.push_back(-1);
    for (size_t i=0; i<remove_catalog.size(); i++) new_conflicts.push_back(i);
      //conflict_group[-1].push_back(i);
    conflict_group.push_back(new_conflicts);
  }
  // return conflict_group;
  // std::cout<<"got conflict_group"<<"\n";
  if(MHS_debug) std::cout<<"created ConflictGroup of size "+Str((int)conflict_group.size())<<endl;;
  if(MHS_debug) std::cout<<PrintConflictGroup(conflict_group,remove_catalog,add_catalog)+"\n";
}
/**

TODO: energy evaluation
      reverse addition deletion of pairs in pair table (only those, which were done (i.e. don't reverse removed pairs
      that were not treated because they already have been dealt with before)

*/



/**
   idx - idx of the conflict group that generates the trail
*/



void PartialPath::UndoPartialPathTrail(){
  if(MHS_debug) cout<<"UndoPartialPathTrail of size "+Str((int)currentPartialPathTrail.size())<<endl;
  // for(size_t i=currentPartialPathTrail.size()-1;i>=0;i--){
    for(std::vector<int>::reverse_iterator it=currentPartialPathTrail.rbegin();it!=currentPartialPathTrail.rend();it++){
      //      int idx=*(it-1);//currentPartialPathTrail[i];
      int idx=*it;
    if(MHS_debug) cout<<"idx "+Str(idx)<<std::endl;
    //base pair was added
    if(idx>=BP_ADD_CONST-1) {
      idx-=BP_ADD_CONST;
      //ignore -1, as it means nothing got added and the index is only for the emptying of the conflict_group
      if(idx==-1) continue;
      RemoveFromPairTable(add_catalog[idx]);
      }
    //base pair was removed
    else{
      AddToPairTable(remove_catalog[idx]);
    }
  }
  currentPartialPathTrail.clear();
}


bool PartialPath::MinusculeStackIsInTarget(Stack mStack){
  if(MHS_debug) cout<<"MinusculeStackIsInTarget with stack of size "+Str((int)mStack.bp.size())<<std::endl;
  for(size_t i=0;i<mStack.bp.size();i++){
    if(MHS_debug) cout<<i+" "+PrintBasePair(mStack.bp[i])<<endl;
    if(find(target.begin(),target.end(),mStack.bp[i])==target.end()) return false;
  }
  return true;

}

int PartialPath::StacksCreatedByBPRemoval(int first,int second){
  
  int count=0;
   if(IsLonelyPair(first,second)) {
     minusculeStack.Clear();
     return 0;
   }
   //the bp is part of the minuscule Stack of size 0. Its removal decreases the size
   //of the minuscule Stack, but doesnt create a new one.
   else if(minusculeStack.GetSize()>0 && minusculeStack.Contains(first,second)){
     return 0;
   }
 //now it can only have been part of a stack that was not minuscule but might become minuscule due
 //to the removal
 //check elts in pairtable around bp to obtain whether the remaining stack is minuscule
 int inside_count=0;
 int inside_first=first+1;
 int inside_second=second-1;
 bool inside_GC=true;
 while(PairIsInPairtable(inside_first,inside_second)){
   if(!Stack::IsGCPair(sequence[inside_first-1],sequence[inside_second-1])) inside_GC=false;
   inside_first++;
   inside_second--;
   inside_count++;
 }
 bool inside_minuscule=Stack::IsMinuscule(inside_count,inside_GC);
 if(inside_count>0 && inside_minuscule) count++;


 //=================================================
 int outside_count=0;
 int outside_first=first-1;
 int outside_second=second+1;
 bool outside_GC=true;
 while(PairIsInPairtable(outside_first,outside_second)){
   if(!Stack::IsGCPair(sequence[outside_first-1],sequence[outside_second-1])) outside_GC=false;
   outside_first--;
   outside_second++;
   outside_count++;
 }
 bool outside_minuscule=Stack::IsMinuscule(outside_count,outside_GC);
 if(outside_count>0 && outside_minuscule) count++;
 //only one minuscule stack
 createdMinusculeStack1.Clear();
 createdMinusculeStack2.Clear();
 if(count==1){ 
   if(inside_minuscule) {
     inside_first--;
     inside_second++;
     if(MHS_debug) std::cout<<"create inside_minuscule Stack with first "+Str(inside_first)<<endl;
     while(inside_first!=first){
       createdMinusculeStack1.AddOutside(inside_first,inside_second,sequence[inside_first-1],sequence[inside_second-1]);
       inside_first--;
       inside_second++;
     }
   }
   else if(outside_minuscule){
     outside_first++;
     outside_second--;
     if(MHS_debug) std::cout<<"create outside_minuscule Stack with first "+Str(outside_first)<<endl;
     while(outside_first!=first){
       createdMinusculeStack1.AddInside(outside_first,outside_second,sequence[outside_first-1],sequence[outside_second-1]);
       outside_first++;
       outside_second--;
     }
   }
   if(MHS_debug) std::cout<<"createdMinusculeStack1 "+createdMinusculeStack1.Print()<<endl;
   if(MinusculeStackIsInTarget(createdMinusculeStack1)){
     if(MHS_debug)  std::cout<<"ignore bc in target "<<endl;
     createdMinusculeStack1.Clear();
     count--;
   }
 }
 if(count==2){
   //add inside to stack1, outside to stack2
  
 
     inside_first--;
     inside_second++;
     while(inside_first!=first){
       createdMinusculeStack1.AddOutside(inside_first,inside_second,sequence[inside_first-1],sequence[inside_second-1]);
       inside_first--;
       inside_second++;
     }
   
   
     outside_first++;
     outside_second--;
     while(outside_first!=first){
       createdMinusculeStack2.AddInside(outside_first,outside_second,sequence[outside_first-1],sequence[outside_second-1]);
       outside_first++;
       outside_second--;
     }

   if(MHS_debug) std::cout<<"createdMinusculeStack1 "+createdMinusculeStack1.Print()<<endl;
   if(MinusculeStackIsInTarget(createdMinusculeStack1)){
      if(MHS_debug) std::cout<<"ignore bc in target "<<endl;
      createdMinusculeStack1.Clear();
      count--;
   }

   if(MHS_debug) std::cout<<"createdMinusculeStack2 "+createdMinusculeStack2.Print()<<endl;
   if(MinusculeStackIsInTarget(createdMinusculeStack2)){
      if(MHS_debug) std::cout<<"ignore bc in target "<<endl;
      createdMinusculeStack2.Clear();
      count--;
   }

 }
 

 return count;

}



/*
bool PartialPath::BPRemovalCreatesSecondMinusculeStack(int first,int second){
  //do not allow creating 2 minuscule stacks at once.
  //  if(minusculeStack.GetSize()==0) return false;
  //if the pair was single, it constituted the minusculeStack, which therefore disappears
 if(IsLonelyPair(first,second)) {
   minusculeStack.Clear();
   return false;
 }
 //now it can only have been part of a stack that was not minuscule but might become minuscule due
 //to the removal
 //check elts in pairtable around bp to obtain whether the remaining stack is minuscule
 int inside_count=0;
 int inside_first=first+1;
 int inside_second=second-1;
 bool inside_GC=true;
 while(PairIsInPairtable(inside_first,inside_second)){
   if(!Stack.IsGCPair(sequence[inside_first],sequence[inside_second])) inside_GC=false;
   inside_first++;
   inside_second--;
   inside_count++;
 }
 if(inside_count>0 && Stack.IsMinuscule(inside_count,inside_GC))
 int outside_count=0;
  return true;
}
*/
void PartialPath::RemoveMinusculeStack(){
  RemoveMinusculeStack(minusculeStack);
  minusculeStack.Clear();
}

void PartialPath::RemoveMinusculeStack(Stack mStack){
  if(MHS_debug) cout<<"RemoveMinusculeStack "+mStack.Print()<<endl;
 
  //record movements in a trail, so they can be undone


  //if the Stack is increasing, add further base pairs to it
  //if tgt contains stack
  //the stack is here to stay if both its inner and outermost bp are in the tgt, otherwise it is not!!!
  if(find(target.begin(),target.end(),mStack.bp.back())!=target.end() &&find(target.begin(),target.end(),mStack.bp.front())!=target.end() ){
  
      //do we add base pairs inwards or outwards? 
      //inwards
      int inwards_first=mStack.bp[0].first+1;
      int inwards_second=mStack.bp[0].second-1;
      //      while(mStack.CanBeAddedInside(inwards_first,inwards_second)) {
      int addCatalogIdx;
      while((addCatalogIdx=find(add_catalog.begin(),add_catalog.end(),std::make_pair(inwards_first,inwards_second))-add_catalog.begin())<(int)add_catalog.size()){
       mStack.AddInside(inwards_first,inwards_second,sequence[inwards_first-1],sequence[inwards_second-1]);
       AddToPairTableWithLog(inwards_first,inwards_second,addCatalogIdx);
       inwards_first++;
       inwards_second--;
       if(!mStack.IsMinuscule()) break;
      }
      //outwards
      int outwards_first=mStack.bp.back().first-1;
      int outwards_second=mStack.bp.back().second+1;
      //while(mStack.CanBeAddedOutside(outwards_first,outwards_second)) {
      while((addCatalogIdx=find(add_catalog.begin(),add_catalog.end(),std::make_pair(outwards_first,outwards_second))-add_catalog.begin())<(int)add_catalog.size()){
       mStack.AddOutside(outwards_first,outwards_second,sequence[outwards_first-1],sequence[outwards_second-1]);
       AddToPairTableWithLog(outwards_first,outwards_second,addCatalogIdx);
       outwards_first--;
       outwards_second++;
       if(!mStack.IsMinuscule()) break;
      }

  }
 
  //if it is decreasing, remove all base pairs from it 
  //if src contains stack
  //  while(mStack.GetSize()>0){
 else{
  for(int i=0;i<mStack.GetSize();i++){
    //do not remove basepairs that are in the target
    if(find(target.begin(),target.end(),mStack.bp[i])!=target.end()) continue;    

    int first=mStack.bp[i].first;
    int second=mStack.bp[i].second;
    if(!PairIsInPairtable(first,second)) continue;

    int removeCatalogIdx=find(remove_catalog.begin(),remove_catalog.end(),mStack.bp[i])-remove_catalog.begin();
    RemoveFromPairTableWithLog(first,second,removeCatalogIdx);
  }
 }

}


bool PartialPath::BPAdditionCreatesSecondMinusculeStack(int first,int second){
 if(MHS_debug) std::cout<<"BPAdditionCreatesSecondMinusculeStack "+Str(first)+" "+Str(second)<<endl;


  if(minusculeStack.GetSize()==0) return false;
  //if we are here, there already is a minuscule Stack
  //the only way to create a new  minuscule Stack is if the bp is lonely
  //and if is still part of mStack in the target structure. As only base pairs in the
  //target structure ever get added, this should never occur.

  //determines whether the stack around a pair (first,second) is minuscule in target.
  if(IsLonelyPair(first,second)) {
    if(MHS_debug) std::cout<<"added pair is lonely"<<endl;
    bool GC=Stack::IsGCPair(sequence[first-1],sequence[second-1]);
    //is in target automatically
    std::pair<int,int> bp=std::make_pair(first,second);
    int stack_size=1;
    //go inside
    std::pair<int,int> bp_inner= std::pair<int,int>(bp.first+1,bp.second-1);
    while(true){
      //finds base pair outside of bp
      if(find(target.begin(),target.end(),bp_inner)!=target.end()){
        stack_size++;
        if(!IsGCPair(bp_inner,sequence)) GC=false;
        bp_inner= std::pair<int,int>(bp_inner.first+1,bp_inner.second-1);
      }
      else break;
    }

    std::pair<int,int> bp_outer= std::pair<int,int>(bp.first-1,bp.second+1);
    while(true){
      //finds base pair outside of bp
      if(find(target.begin(),target.end(),bp_outer)!=target.end()){
        stack_size++;
        if(!IsGCPair(bp_outer,sequence)) GC=false;
        bp_outer= std::pair<int,int>(bp_outer.first-1,bp_outer.second+1);
      }
      else break;
    }
    //if the stack is minuscule in the target structure, it is not counted as minuscule
    //contribution to the current structure, so we return false
    if(Stack::IsMinuscule(stack_size,GC)) return false;
    else{
      if(MHS_debug) std::cout<<"and creates new mStack"<<endl;
      return true;
    }
  }
  //the minusculeStack disappears if the current bp increases its size so it is not minuscule any more
  
  if(minusculeStack.CanBeAddedOutside(first,second)) 
    minusculeStack.AddOutside(first,second,sequence[first-1],sequence[second-1]);
  else if(minusculeStack.CanBeAddedInside(first,second)) 
    minusculeStack.AddInside(first,second,sequence[first-1],sequence[second-1]);
  if(!minusculeStack.IsMinuscule()) minusculeStack.Clear();

  //the count of minuscule Stacks is now either 0 or 1, so we didnt create a second one.
  return false;
}



bool PartialPath::IsLonelyPair(int first,int second){
  bool no_bp_outside= (pair_table[first-1]!=second+1);
  bool no_bp_inside= (pair_table[first+1]!=second-1);
  if(no_bp_outside && no_bp_inside) return true;
  else return false;
}


void PartialPath::TreatConflictGroupIndex(int idx){
  bool debug=false;
 if(MHS_debug) cout<<"TreatConflictGroupIndex"<<endl;
  if(debug) cout<<PrintPairTable()<<endl;
 //the items to remove are stored from the second element onwards.
  if(MHS_debug) cout<<"print conflict_group["+Str(idx)+"]"<<endl;
  //for(size_t i=1;i<conflict_group[idx].size();i++){
   
  // }
  for(size_t i=1;i<conflict_group[idx].size();i++){
    //    int remove_idx=conflict_group[idx][i];
   
    //if(!RemovePairIsPresent(idx,remove_idx)) continue;
    int first=GetFirstFromConflictGroupRemove(idx,i);
    int second=GetSecondFromConflictGroupRemove(idx,i);
    // int first=remove_catalog[remove_idx].first;//GetFirstFromConflictGroupRemove(idx,remove_idx);
    //    int second=remove_catalog[remove_idx].second;//GetSecondFromConflictGroupRemove(idx,remove_idx);
    if(MHS_debug) cout<<"remove idx "+Str(conflict_group[idx][i])+" bp "+PrintBasePair(std::make_pair(first,second))<<endl;
    //do nothing if the pair has already been removed from pair_table
    if(!PairIsInPairtable(first,second)) {
        if(MHS_debug)cout<<"Is not in pair_table"<<endl;     
        continue;
    }
    if(MHS_debug) {
    if(minusculeStack.GetSize()>0) cout<<"Already have a minuscule stack"<<endl;
    else cout<<"Dont have a minuscule stack"<<endl;
    }
    int createdmStackCount=StacksCreatedByBPRemoval(first,second);
    int totalmStackCount=createdmStackCount;
   
    if(minusculeStack.GetSize()>0) totalmStackCount++;
    if(MHS_debug) cout<<"createdmStackCount "+Str(createdmStackCount)+" totalmStackCount "+Str(totalmStackCount)<<endl;
    if(totalmStackCount==1 && createdmStackCount==1){
      if(createdMinusculeStack1.GetSize()>0) minusculeStack=createdMinusculeStack1;
      else if(createdMinusculeStack2.GetSize()>0) minusculeStack=createdMinusculeStack2;
    }


    if(totalmStackCount>=2){
      //if there already was a mStack, remove it
      if(totalmStackCount>createdmStackCount){
        if(MHS_debug) cout<<"just created a mStack, so remove the previous mStack"<<endl;
        RemoveMinusculeStack();
       //determine the new minuscule stack
      }
      //if two minusculeStacks just got created, remove the smaller of the two and keep the larger of the two as minusculeStack 
     //remove the smaller of the two created mStacks, and make the bigger the new mStack
     if(createdmStackCount>=1) {
       if(MHS_debug) cout<<" createdMinusculeStack1.GetSize() "+Str(createdMinusculeStack1.GetSize())+" createdMinusculeStack2.GetSize() "+Str(createdMinusculeStack2.GetSize())<<endl;
       if(createdMinusculeStack1.GetSize()>createdMinusculeStack2.GetSize()){
          if(createdMinusculeStack2.GetSize()>0) RemoveMinusculeStack(createdMinusculeStack2);
          minusculeStack=createdMinusculeStack1;
       }
       else if(createdMinusculeStack2.GetSize()>createdMinusculeStack1.GetSize()){
          if(createdMinusculeStack1.GetSize()>0) RemoveMinusculeStack(createdMinusculeStack1);
          minusculeStack=createdMinusculeStack2;
       }
      if(MHS_debug) std::cout<<"MinusculeStack after removal "+minusculeStack.Print()<<endl;
     }
    }
    int removeCatalogIdx=GetRemoveCatalogIndex(idx,i);//remove_idx);
    RemoveFromPairTableWithLog(first,second,removeCatalogIdx);

  }
  //if to_add is -1, all has already been added, so what comes below is skipped.
  if(conflict_group[idx][0]==-1) {
    //it needs to be recorded in the trail, so the conflict_group can be emptied.
    currentPartialPathTrail.push_back(-1+BP_ADD_CONST);
    return;
  }
  //the element at given index inthe conflictgroup
  int first=GetFirstFromConflictGroupAdd(idx);
  int second=GetSecondFromConflictGroupAdd(idx);
  int addCatalogIdx=GetAddCatalogIndex(idx);
  if(BPAdditionCreatesSecondMinusculeStack(first,second)) {
     RemoveMinusculeStack();
     minusculeStack.AddOutside(first,second,sequence[first-1],sequence[second-1]);
  }
  //add the base pair to pair_table unless it is already part of it;
  if(!PairIsInPairtable(first,second)) AddToPairTableWithLog(first,second,addCatalogIdx);

  }


/*
    std::vector<Stack> newMinusculeStacks=GetMinusculeStacksForRemoval(first,second);

       if((int)newMinusculeStacks.size()>=minusculeMaxCount){
       //1.A second stack jsut appeared via pushback
       //2.remove the original stack, i.e. the first one in newMinusculeStacks
       Stack firstMinuscule=minusculeStacks[0];
       //get rid of all elements in the minuscule Stack.
       //if its elements will be part of the final structure, we have to add bp until the stack is larger than minuscule
       //if they won't be part, we have to remove all of its elements.
       if(firstMinuscule.persists){

       }
       else{


       }
       minusculeStacks.erase(minusculeStacks.begin());
       }
*/



bool PartialPath::PairIsInPairtable(int first,int second){
    if(first<=0) return false;
    bool firstPaired=(pair_table[first]==second);
    bool secondPaired=(pair_table[second]==first);
    return (firstPaired && secondPaired);
 }



/*
 bool PartialPath::RemovePairIsPresent(int first,int second){
    bool firstPaired=(pair_table[first]==second);
    bool secondPaired=(pair_table[second]==first);
    return (firstPaired && secondPaired);
 }
 */

 /*
 bool PartialPath::RemovePairIsPresent(int group_idx,int remove_idx){
    int first=GetFirstFromConflictGroupRemove(group_idx,remove_idx);
    int second=GetSecondFromConflictGroupRemove(group_idx,remove_idx);
    bool firstPaired=(pair_table[first]==second);
    bool secondPaired=(pair_table[second]==first);
    return (firstPaired && secondPaired);
  }
 */

 

  int PartialPath::GetFirstFromConflictGroupRemove(int group_idx,int remove_idx){
    bool debug=false;
    if(debug){
      cout<<"GetFirstFromConflictGroupRemove group member "+Str(group_idx)+" idx in list "+Str(remove_idx)<<endl;
      cout<<"yields index in remove catalog "+Str(conflict_group[group_idx][remove_idx])<<endl;
      cout<<"bp "+PrintBasePair(remove_catalog[conflict_group[group_idx][remove_idx]])<<endl;
    }
    return remove_catalog[conflict_group[group_idx][remove_idx]].first;
  } 

  int PartialPath::GetSecondFromConflictGroupRemove(int group_idx,int remove_idx){
    bool debug=false;
    if(debug){
      cout<<"GetSecondFromConflictGroupRemove group member "+Str(group_idx)+" idx in list "+Str(remove_idx)<<endl;
      cout<<"yields index in remove catalog "+Str(conflict_group[group_idx][remove_idx])<<endl;
      cout<<"bp "+PrintBasePair(remove_catalog[conflict_group[group_idx][remove_idx]])<<endl;
    }
    return remove_catalog[conflict_group[group_idx][remove_idx]].second;
  } 


  int PartialPath::GetAddCatalogIndex(int idx){
    return conflict_group[idx][0];
  }



  int PartialPath::GetRemoveCatalogIndex(int group_idx,int remove_idx){
    if(MHS_debug) cout<<"GetRemoveCatalogIndex group member "+Str(group_idx)+" idx in list "+Str(remove_idx)+" ";
    if(MHS_debug) cout<<"yields index in remove catalog "+Str(conflict_group[group_idx][remove_idx])<<endl;
    return conflict_group[group_idx][remove_idx];
  }




 int PartialPath::GetFirstFromConflictGroupAdd(int group_idx){
    return add_catalog[conflict_group[group_idx][0]].first;
  } 

  int PartialPath::GetSecondFromConflictGroupAdd(int group_idx){
    return add_catalog[conflict_group[group_idx][0]].second;
  } 

  void PartialPath::AddToPairTable(std::pair<int,int> bp){
    // cout<<"AddToPairTable"+PrintBasePair(bp) <<std::endl;
    pair_table[bp.first]=bp.second;
    pair_table[bp.second]=bp.first;
  }


  void PartialPath::RemoveFromPairTable(std::pair<int,int> bp){
    // cout<<"RemoveFromPairTable" +PrintBasePair(bp)<<std::endl;
    pair_table[bp.first]=0;
    pair_table[bp.second]=0;
  }


  void PartialPath::AddToPairTable(int first,int second){
    //    cout<<"AddToPairTable" +PrintBasePair(std::make_pair(first,second))<<std::endl;
    pair_table[first]=second;
    pair_table[second]=first;
  }



  void PartialPath::RemoveFromPairTable(int first,int second){
    //   cout<<"RemoveFromPairTable"+PrintBasePair(std::make_pair(first,second)) <<std::endl;
    pair_table[first]=0;
    pair_table[second]=0;
  }

  void  PartialPath::RemoveFromPairTableWithLog(int first,int second,int removeCatalogIdx){
    //    cout<<"RemoveFromPairTableWithLog idx "+Str(removeCatalogIdx)<<endl;
    currentPartialPathTrail.push_back(removeCatalogIdx);
    RemoveFromPairTable(first,second);
  }

 void  PartialPath::AddToPairTableWithLog(int first,int second,int addCatalogIdx){
   //   cout<<"AddToPairTableWithLog idx "+Str(addCatalogIdx)<<endl;
    currentPartialPathTrail.push_back(addCatalogIdx+BP_ADD_CONST);
    AddToPairTable(first,second);
  }



PartialPath::PartialPath(std::vector<std::pair<int,int> > src,std::vector<std::pair<int,int> > target,std::string sequence,int lookahead,const std::vector<std::pair<int,int> > & remove_catalog,const std::vector<std::pair<int,int> > & add_catalog){
  this->src=src;
  this->target=target;
  this->sequence=sequence;
  this->lookahead=lookahead;
  this->add_catalog =add_catalog;
  this->remove_catalog =remove_catalog;
}


double PartialPath::WalkPartialPath( std::vector<int> combination ){

 if(MHS_debug)  cout<<"WalkPartialPath \n";
  // std::string current_structure=sequence; 
  //std::string highest_structure;

  //  double highest_energy =  std::numeric_limits<double>::min();
    double highest_energy =  -INF;//std::numeric_limits<double>::min();
  //double current_value=std::numeric_limits<double>::min();  
  



  for (size_t k=0; k<combination.size(); k++) { 
    int idx = combination[k]-1;
    TreatConflictGroupIndex(idx);
    double current_energy=FastEvalEnergy(sequence);
    if(MHS_debug) cout<<"current_energy "+Str(current_energy)<<endl;
    if(current_energy>highest_energy) highest_energy=current_energy;      
    // UndoPartialPathTrail();
  }
 if(MHS_debug)  cout<<"currentPartialPathTrail after WalkPartialPath "+PrintCombination(currentPartialPathTrail)+"\n";
    return highest_energy;
}
    /**
    std::map<int,std::vector<int> >::const_iterator it=conflict_group.begin();
    advance(it,idx);
    std::pair<int,std::vector<int> > conflict=*it;
    //XXXXXXXXXXXXXXXXXXXXXXXXX DO WE NEED TO TREAT TO AD HERE?
    //    int to_add=conflict.first;
    //if(to_add==1) continue;
    std::vector<int> to_remove=conflict.second;
    // remove what impedes adding the pair in combination[k], then
    // add it.  TREAT ONE ELEMENT IN THE CONFLICT GROUP
    for (size_t l=0; l<to_remove.size(); l++) {

     std::pair<int,int> bp_to_remove=only_in_base_pairs[to_remove[l]]
     std::vector<Stack> newMinusculeStacks= PartialPath::GetMinusculeStacksForRemoval(bp_to_remove);
     if(newMinusculeStacks.size()>=minusculeMaxCount){
       //1.A second stack jsut appeared via pushback
       //2.remove the original stack, i.e. the first one in newMinusculeStacks
       Stack firstMinuscule=minusculeStacks.pop_front();//[0];
       //get rid of all elements in the minuscule Stack.
       //if its elements will be part of the final structure, we have to add bp until the stack is larger than minuscule
       //if they won't be part, we have to remove all of its elements.
       if(firstMinuscule.persists){

       }
       else{


       }
}
    

    





      //===========================================================
      //std::vector<std::pair<int,int> >::iterator it = remove(current.begin(),current.end(),only_in_base_pairs[to_remove[l]]);
      //if only_in_base_pairs[to_remove[l] hasn't been removed before
      if(it!=current.end()){
         std::vector<std::pair<int,int> > current_bp_after_remove= std::vector<std::pair<int,int> >(current.begin(),it);

        std::vector<std::pair<int,int> > lonely_bp=GetMinusculeStacks();
	//        std::vector<std::pair<int,int> > lonely_bp=GetMinusculeStacks(current_bp_after_remove,base_pairs2,sequence);

        int lonely_bp_count=lonely_bp.size();
        //if removing the requested base pair results in too many lonely base pairs, we have to get rid of the initial lonely base pair first
        if(lonely_bp_count>minusculeMaxCount){
          //need to keep the original lonely_bp in memeory
	}

         pair_table[only_in_base_pairs[to_remove[l]].first]=0;
         pair_table[only_in_base_pairs[to_remove[l]].second]=0;
         removed_pairs.push_back(to_remove[l]);
      }
     
      current.erase(it, current.end());
      //===========================================================


    }
    */
    //}
    
  


/*
void PartialPath::RemoveFromConflictGroup(std::vector<int> added_regular,std::vector<std::pair<int,int> > added_correctives){

  for(size_t i=0;i<added_regular.size();i++){
    conflict_group.erase(added_regular[i]);
  }
  //bp that are added are in  only_in_base_pairs2, so this is the container that is searched for the elements of added_correctives;
  for(size_t i=0;i<added_correctives.size();i++){
    std::vector<std::pair<int,int> >::iterator added_corrective=find(only_in_base_pairs2.begin(),only_in_base_pairs2.end(),added_correctives[i]);
    if(added_corrective!=only_in_base_pairs2.end()){
      int idx=added_corrective-only_in_base_pairs2.begin();
      conflict_group.erase(idx);
    }
  }

}
*/
/*
// CONSIDER THE POSSIBLE COMBINATIONS OF THE GIVEN CONFLICT GROUP IN TURN. 
  for (size_t k=0; k<combination.size(); k++) { 
    int idx = combination[k]-1;
    std::map<int,std::vector<int> >::const_iterator it=conflict_group.begin();
    advance(it,idx);
    std::pair<int,std::vector<int> > conflict=*it;
    int to_add=conflict.first;
    //if(to_add==1) continue;
    std::vector<int> to_remove=conflict.second;
    // remove what impedes adding the pair in combination[k], then
    // add it.  TREAT ONE ELEMENT IN THE CONFLICT GROUP
    for (size_t l=0; l<to_remove.size(); l++) {
      size_t old_size=current_pairs.size();
      // TREAT ONE BASE PAIR
      #ifdef _DEBUG_MH_
	  std::cout<<"remove pair ("+Str(only_in_base_pairs[to_remove[l]].first)+","+Str(only_in_base_pairs[to_remove[l]].second)+")"<<std::endl;
      #endif
      std::vector<std::pair<int,int> >::iterator it = remove(current_pairs.begin(),current_pairs.end(),only_in_base_pairs[to_remove[l]]);


      //if only_in_base_pairs[to_remove[l] hasn't been removed before
      if(it!=current_pairs.end()){
         std::vector<std::pair<int,int> > current_bp_after_remove= std::vector<std::pair<int,int> >(current_pairs.begin(),it);

        std::vector<std::pair<int,int> > lonely_bp=GetMinusculeStacks(current_bp_after_remove,base_pairs2,sequence);
        int lonely_bp_count=lonely_bp.size();
        //if removing the requested base pair results in too many lonely base pairs, we have to get rid of the initial lonely base pair first
        if(lonely_bp_count>max_lonely_bp){
          //need to keep the original lonely_bp in memeory
	}

         pair_table[only_in_base_pairs[to_remove[l]].first]=0;
         pair_table[only_in_base_pairs[to_remove[l]].second]=0;
         removed_pairs.push_back(to_remove[l]);
      }
      current_pairs.erase(it, current_pairs.end());
*/
void PartialPath::AddBestPartialPathTrailToPathTrail(){
 if(MHS_debug)  cout<<"AddBestPartialPathTrailToPathTrail \n";

  for(size_t i=0;i<bestPartialPathTrail.size();i++){
    int idx=bestPartialPathTrail[i];
    pathTrail.push_back(idx);
    if(idx>=BP_ADD_CONST-1){
      idx-=BP_ADD_CONST;
      //this is new
      if(idx==-1) continue;
      //add the bp to the pair_table and latest structure
      std::pair<int,int> bp=add_catalog[idx];
      pair_table[bp.first]=bp.second;
      pair_table[bp.second]=bp.first;
    }
    else{
      //remove the bp from the pair_table and latest structure
      std::pair<int,int> bp=remove_catalog[idx];
      pair_table[bp.first]=0;
      pair_table[bp.second]=0;
    }
  
  }
  //  bestPartialPathTrail.clear();
if(MHS_debug)   cout<<"yields pathTrail "+ PrintCombination(pathTrail)+"\n";
}


void PartialPath::FindOptimalPath(std::vector<std::pair<double,std::string> >& structures){
  bool debug=false;
  if(MHS_debug)cout<<"FindOptimalPath\n";



 
  std::vector<int> combination= std::vector<int>();  

   
  GetConflictGroup();
  
  while(conflict_group.size()>0){
    double best_combination_saddle_energy = INF;//std::numeric_limits<double>::max();
   long int n_combinations = N_take_k( conflict_group.size(),lookahead);


  std::vector<Stack> minusculeStacks=GetMinusculeStacks();
  //if there is more than one minuscule Stack, rremove all but the largest.
  int greatestMStackIndex=-1;
  int maxsize=0;  
  for(size_t i=0;i<minusculeStacks.size();i++){
    int size=minusculeStacks[i].GetSize();
    if(size>maxsize) {
       maxsize=size;
       greatestMStackIndex=i;
    }
  }
  Stack originalMinusculeStack;
  if(greatestMStackIndex==-1) originalMinusculeStack=Stack();
  else                        originalMinusculeStack=minusculeStacks[greatestMStackIndex];
  //get rid of the smaller stacks
  for(size_t i=0;i<minusculeStacks.size();i++){
    if((int)i==greatestMStackIndex ) continue;
    RemoveMinusculeStack(minusculeStacks[i]);
  }


   //  cout<<"Treat conflict group of size"+Str((int)conflict_group.size())+"\n";
    //try all combinations and record combination idx with the best energy.Keep track of saddle point for the entire function all, that is throughout all groups.
    for (long int j=0; j<n_combinations; j++) {
      minusculeStack=originalMinusculeStack;

      if(MHS_debug)cout<<"combination "+Str((int)j)+" out of "+Str((int)n_combinations)<<endl;
      combination = GetCombination(conflict_group.size(),lookahead,j);
      double saddle_energy=WalkPartialPath(combination);
        if(MHS_debug)cout<<"obtains saddle "+Str(saddle_energy)+"\n";
      if(saddle_energy<best_combination_saddle_energy) {
        if(debug) cout<<"UpdateBestPartialPathTrail from current Trail of size "+Str((int)currentPartialPathTrail.size())+"\n";
        best_combination_saddle_energy=saddle_energy;
        bestPartialPathTrail=currentPartialPathTrail;
      }
      UndoPartialPathTrail();
    }
    if(debug) cout<<"AddBestPartialPathTrail of size "+Str((int)bestPartialPathTrail.size())+" ToPathTrail"<<std::endl;
    //updates pair_table
    AddBestPartialPathTrailToPathTrail();
    //remove the used up elements from the conflict group
    if(MHS_debug)cout<<"remove the used up elements from the conflict group of size "+Str((int)conflict_group.size())<<std::endl;
    if(MHS_debug) std::cout<<PrintConflictGroup(conflict_group,remove_catalog,add_catalog)<<std::endl;

    for(size_t i=0;i<bestPartialPathTrail.size();i++){
     if(debug) cout<<"bestPartialPathTrail ["+Str((int)i)+"] "+Str(bestPartialPathTrail[i])<<std::endl;
     if(bestPartialPathTrail[i]<BP_ADD_CONST-1) continue;
     
     //do not change horses midstream: ain't happening, as we break as soon as an element is deleted
     for(size_t j=0;j<conflict_group.size();j++){
       if(conflict_group[j][0]==bestPartialPathTrail[i]-BP_ADD_CONST){
         conflict_group.erase(conflict_group.begin()+j);
	 if(debug)std::cout<<"erase "+Str((int)j)<<endl;
         break;
       }
     }
     if(debug)cout<<"conflict_group has now size "+Str((int)conflict_group.size())<<std::endl;
    }

  }
  //call this before UpdateBacktrackBase, while add_catalog and remove_catalog are not emptied yet.
  if(debug) cout<<"AddPathTrailToStructureList"<<std::endl;
  AddPathTrailToStructureList(structures);
  UpdateBacktrackBase();

  //cout<<std::endl;
  pathTrail.clear();
}

void PartialPath::AddPathTrailToStructureList(std::vector<std::pair<double,std::string> >& structures){
  bool debug=false;
if(debug)  cout<<PrintPairTable()<<endl;
  std::string latest_structure=structures.back().second;
  MakePairTable(const_cast<char*>(latest_structure.c_str()));
  //turn latest structure into pair_table
  for(size_t i=0;i<pathTrail.size();i++){
    int idx=pathTrail[i];
   if(debug) cout<<"idx "<<idx<<endl;
    if(idx>=BP_ADD_CONST-1){
      idx-=BP_ADD_CONST;
      //-1 indicates nothing was added, it only holds the place in the conflict group
      if(idx==-1) continue;
      //add the bp to the pair_table and latest structure
      std::pair<int,int> bp=add_catalog[idx];
      latest_structure[bp.first-1]='(';
      latest_structure[bp.second-1]=')';
      pair_table[bp.first]=bp.second;
      pair_table[bp.second]=bp.first;
      double energy=FastEvalEnergy(sequence);
      structures.push_back(std::make_pair(energy,latest_structure));
    }
    else{
      //remove the bp from the pair_table and latest structure
      std::pair<int,int> bp=remove_catalog[idx];
      latest_structure[bp.first-1]='.';
      latest_structure[bp.second-1]='.';
      pair_table[bp.first]=0;
      pair_table[bp.second]=0;
      double energy=FastEvalEnergy(sequence);
      structures.push_back(std::make_pair(energy,latest_structure));

    }
     if(debug) cout<<"becomes idx "<<idx<<endl;
  }
  if(debug){
    cout<<"New Structure List:"<<endl;
    for(size_t l=0;l<structures.size() ;l++){
      cout<<structures[l].second+" "+Str(structures[l].first)<<endl;
    }
    cout<<PrintPairTable()<<endl;
  }
  //  pathTrail.clear();
}



//erase every index listed in pathTrail from the add/remove catalog
void PartialPath::UpdateBacktrackBase(){
  bool debug=false;
       if(MHS_debug) cout<<"UpdateBacktrackBase with pathTrail of size "+Str((int)pathTrail.size())<<std::endl;
       
       for (int i=0; i<(int)pathTrail.size(); i++) {
          if(debug) cout<<"pathTrail["+Str(i)+"] "+Str(pathTrail[i])<<endl;
      }
       if(MHS_debug) cout<<"add_catalog size "+Str((int)add_catalog.size())+" remove_catalog size "+Str((int)remove_catalog.size())<<std::endl;
       //the lowest index comes at the end
      sort(pathTrail.begin(),pathTrail.end(),std::greater<int>());
      for (int i=0; i<(int)pathTrail.size(); i++) { 
	int idx = pathTrail[i];
if(MHS_debug) 	cout<<"idx "+Str(idx)+" ";
        if(idx>=BP_ADD_CONST-1){
	  idx-=BP_ADD_CONST;
          //-1 is place holder in conflict group but does not affect add_catalog or src, as it is not associated with a base pair
          if(idx==-1) continue;
           std::vector<std::pair<int,int> >::iterator it=add_catalog.begin();
          advance(it,idx);
     if(MHS_debug)      cout<<"erase "+PrintBasePair(*it)<<endl;
          add_catalog.erase(it);
          src.push_back(add_catalog[idx]);

	}
        else {
          std::vector<std::pair<int,int> >::iterator it=remove_catalog.begin();
          advance(it,idx);
       if(MHS_debug)    cout<<"erase "+PrintBasePair(*it)<<endl;
          remove_catalog.erase(it);
          std::vector<std::pair<int,int> >::iterator it2=src.begin();
          advance(it2,idx);
          src.erase(it2);

	}
      }
  if(MHS_debug)  cout<<"become add_catalog size "+Str((int)add_catalog.size())+" remove_catalog size "+Str((int)remove_catalog.size())<<endl;
}   








  


/*
Returns lonely basepairs in the current conformation that are not also contained in the target confromation
For stacks of length two that are not GC stackings, the outer base pair is returned.
Assumes that base pairs are ordered increasingly in their first digits.
*/


std::vector<Stack> PartialPath::GetMinusculeStacks(){
  bool debug=false;

  if(debug)   std::cout<<"GetMinusculeStacks"<<endl;
  std::vector<Stack> ret=  std::vector<Stack>();
  if(debug) cout<<PrintPairTable()<<endl;
  std::vector< std::pair<int,int> > current=PairTableToBasePairList(pair_table);
  if(debug) std::cout<<"obtained base pair list "+PrintBasePairList(current)<<std::endl;




  for(size_t i=0;i<current.size();i++){
    std::pair<int,int> bp=current[i];
    bool done=false;
    //continue of the current base pair (or a neighboring pair) is already listed in ret 
    //(bc it is the outer bp ina stack, for which the inner bp was already treated)
    for(size_t j=0;j<ret.size();j++){
      if(ret[j].Contains(bp) || ret[j].Contains(bp.first-1,bp.second+1) || ret[j].Contains(bp.first+1,bp.second-1)) {
        done=true;
        break;
      }
    }
    if(done) continue;
    //ignore if this base pair is also part of the final structure
    if(find(target.begin(),target.end(),bp) !=target.end()) continue;
    
    //build stacks by checking whether bases nieghboring bp are paired in the current structure
    bool GC=IsGCPair(bp,sequence);
    int stack_size=1;
    
    std::pair<int,int> bp_outer= std::pair<int,int>(bp.first-1,bp.second+1);
    if(debug) std::cout<<"orig bp_outer "+PrintBasePair(bp_outer)<<endl;
    while(Stack::IsMinuscule(stack_size,GC)  ){
      //finds base pair outside of bp
      if(find(current.begin(),current.end(),bp_outer)!=current.end()){
        stack_size++;
        if(!IsGCPair(bp_outer,sequence)) GC=false;
        bp_outer= std::pair<int,int>(bp_outer.first-1,bp_outer.second+1);
        if(debug) std::cout<<"then bp_outer "+PrintBasePair(bp_outer)<<endl;
      }
      else break;
     
    }

    std::pair<int,int> bp_inner= std::pair<int,int>(bp.first+1,bp.second-1);
    if(debug) std::cout<<"orig bp_inner "+PrintBasePair(bp_inner)<<endl;
    while(Stack::IsMinuscule(stack_size,GC)){
      //finds base pair outside of bp
      if(find(current.begin(),current.end(),bp_inner)!=current.end()){
        stack_size++;
        if(!IsGCPair(bp_inner,sequence)) GC=false;
        bp_inner= std::pair<int,int>(bp_inner.first+1,bp_inner.second-1);
        if(debug) std::cout<<"then bp_inner "+PrintBasePair(bp_inner)<<endl;
      }
      else break;
    }
    if(debug) std::cout<<"stack_size   "<<stack_size<<endl;
    if(!Stack::IsMinuscule(stack_size,GC)) continue;
    if(debug) std::cout<<"is minuscule "<<stack_size<<endl;
    //we actually have a minuscule stack
    Stack mStack=Stack();    
    bp_inner.first--;
    bp_inner.second++;
    if(debug) std::cout<<"bp_inner "+PrintBasePair(bp_inner)<<endl;

     //need to record the outermost pair, which is obtained by moving back in from the last attempt to find an outer
    //base pair, which was too far out.
    // bp_outer= std::pair<int,int>(bp_outer.first+1,bp_outer.second-1);
    //  std::cout<<"bp_outer "+PrintBasePair(bp_outer)<<endl;
    while(bp_inner!=bp_outer){
      mStack.AddOutside(bp_inner.first,bp_inner.second,sequence[bp_inner.first-1],sequence[bp_inner.second-1]);
      //bp_inner now graduall moves OUTSIDE (it moved insde in the loop further up)
      bp_inner= std::pair<int,int>(bp_inner.first-1,bp_inner.second+1);
      if(debug) std::cout<<"new bp_inner "+PrintBasePair(bp_inner)<<endl;

    }
    if(debug) std::cout<<"add stack   "+mStack.Print()<<endl;
    //add the stack if it is minuscule
    ret.push_back(mStack);
  } 
  //std::cout<<"   "<<endl;
  if(debug) std::cout<<"GetMinusculeStacks returns"<<endl;
  return ret;

}




/*
Returns lonely basepairs in the current conformation that are not also contained in the target confromation
For stacks of length two that are not GC stackings, the outer base pair is returned.
Assumes that base pairs are ordered increasingly in their first digits.


std::vector<std::pair<int,int> > PartialPath::GetMinusculeStacks(){

  std::vector<std::pair<int,int> > ret=  std::vector<std::pair<int,int> >();
  for(size_t i=0;i<current.size();i++){
    std::pair<int,int> bp=current[i];
    //continue of the current base pair is already listed in ret (bc it is the outer bp ina stack, for which the inner bp was already treated)
    if(find(ret.begin(),ret.end(),bp)!=ret.end()) continue;
    //ignore if this base pair is also part of the final structure
    if(find(target.begin(),target.end(),bp) !=target.end()) continue;
    
    bool GC=IsGCPair(bp,sequence);
    int stack_size=1;
   
    std::pair<int,int> bp_outer= std::pair<int,int>(bp.first-1,bp.second+1);
    while(true){
      //finds base pair outside of bp
      if(find(current.begin(),current.end(),bp_outer)!=current.end()){
        stack_size++;
        if(!IsGCPair(bp_outer,sequence)) GC=false;
        bp_outer= std::pair<int,int>(bp_outer.first-1,bp_outer.second+1);
      }
      else break;
    }

    std::pair<int,int> bp_inner= std::pair<int,int>(bp.first+1,bp.second-1);
    while(true){
      //finds base pair outside of bp
      if(find(current.begin(),current.end(),bp_inner)!=current.end()){
        stack_size++;
        if(!IsGCPair(bp_inner,sequence)) GC=false;
        bp_inner= std::pair<int,int>(bp_inner.first+1,bp_inner.second-1);
      }
      else break;
    }
    //need to record the outermost pair, which is obtained by moving back in from the last attempt to find an outer
    //base pair, which was too far out.
    bp_outer= std::pair<int,int>(bp_outer.first+1,bp_outer.second-1);

    //if the stack is of size 2, we might already have added it to ret, when the other basepair in the stack was considered
    //in this case we don't add it again.
    if(find(ret.begin(),ret.end(),bp_outer)!=ret.end()) continue;
    //GC stacks need minimal size of 2
    if(GC && stack_size>=2) ret.push_back(bp_outer);
    //other stackings require a size of at least 3
    else if(!GC && stack_size>=3) ret.push_back(bp_outer);


  }
  return ret;

}
*/

std::vector<std::pair<int,int> > base_pairs2=std::vector<std::pair<int,int> > ();

int max_lonely_bp=1;

//std::vector<std::pair<double,std::string> > 
void DoPartialPath(std::vector<std::pair<double,std::string> > & path, const std::vector<int> & combination, std::string sequence, const std::vector<std::pair<int,int> > & 
backtrack_base, const std::map<int,std::vector<int> > & conflict_group, const std::vector<std::pair<int,int> > & only_in_base_pairs, 
const std::vector<std::pair<int,int> > & only_in_base_pairs2){
  /*
void PartialPath(std::vector<std::pair<double,std::string> > & path,const std::vector<int> & combination, std::string sequence,const std::vector<std::pair<int,int> > & backtrack_base, 
const std::map<int,std::vector<int> > & conflict_group,const std::vector<std::pair<int,int> > & only_in_base_pairs,const std::vector<std::pair<int,int> > & only_in_base_pairs2){
  */
  path.clear();
    //  std::vector<std::pair<double,std::string> > path=std::vector<std::pair<double,std::string> > ();

  //  std::cout<<"PartialPath Original PairTable "+PrintPairTable()<<std::endl;
  //bool debug=true;//false;
  std::string highest_structure;
  double highest =  -INF;//std::numeric_limits<double>::min();
  double current_value=-INF;//std::numeric_limits<double>::min();  
  std::string current_structure=sequence;
  std::vector<std::pair<int,int> > current_pairs =backtrack_base;    
  std::vector<int> added_pairs= std::vector<int> ();
  std::vector<int> removed_pairs= std::vector<int> ();
  // CONSIDER THE POSSIBLE COMBINATIONS OF THE GIVEN CONFLICT GROUP IN TURN. 
  for (size_t k=0; k<combination.size(); k++) { 
    int idx = combination[k]-1;
    std::map<int,std::vector<int> >::const_iterator it=conflict_group.begin();
    advance(it,idx);
    std::pair<int,std::vector<int> > conflict=*it;
    int to_add=conflict.first;
    //if(to_add==1) continue;
    std::vector<int> to_remove=conflict.second;
    // remove what impedes adding the pair in combination[k], then
    // add it.  TREAT ONE ELEMENT IN THE CONFLICT GROUP
    for (size_t l=0; l<to_remove.size(); l++) {
      size_t old_size=current_pairs.size();
      // TREAT ONE BASE PAIR
      #ifdef _DEBUG_MH_
	  std::cout<<"remove pair ("+Str(only_in_base_pairs[to_remove[l]].first)+","+Str(only_in_base_pairs[to_remove[l]].second)+")"<<std::endl;
      #endif
      std::vector<std::pair<int,int> >::iterator it = remove(current_pairs.begin(),current_pairs.end(),only_in_base_pairs[to_remove[l]]);


      //if only_in_base_pairs[to_remove[l] hasn't been removed before
      if(it!=current_pairs.end()){
         std::vector<std::pair<int,int> > current_bp_after_remove= std::vector<std::pair<int,int> >(current_pairs.begin(),it);

	 /*
        std::vector<std::pair<int,int> > lonely_bp=GetMinusculeStacks(current_bp_after_remove,base_pairs2,sequence);
        int lonely_bp_count=lonely_bp.size();
        //if removing the requested base pair results in too many lonely base pairs, we have to get rid of the initial lonely base pair first
        if(lonely_bp_count>max_lonely_bp){
          //need to keep the original lonely_bp in memeory
	}
	 */
         pair_table[only_in_base_pairs[to_remove[l]].first]=0;
         pair_table[only_in_base_pairs[to_remove[l]].second]=0;
         removed_pairs.push_back(to_remove[l]);
      }
      current_pairs.erase(it, current_pairs.end());
      if(current_pairs.size()<old_size){
        current_value = FastEvalEnergy(sequence);
	std::string current_structure=  BasePairListToStructure1(sequence.length(),current_pairs);
        path.push_back(std::pair<double,std::string>(current_value,current_structure));
      }
      if ( current_value > highest ) {
	highest = current_value;
	highest_structure = current_structure;
      }
    }
    //-1 means that there is nothing to add any more so only the final contacts are removed
    if(to_add!=-1){
      #ifdef _DEBUG_MH_
          std::cout<<"add pair ("+Str(only_in_base_pairs2[to_add].first)+","+Str(only_in_base_pairs2[to_add].second)+")"<<std::endl;
      #endif
      current_pairs.push_back(only_in_base_pairs2[to_add]);
      pair_table[only_in_base_pairs2[to_add].first]=only_in_base_pairs2[to_add].second;
      pair_table[only_in_base_pairs2[to_add].second]=only_in_base_pairs2[to_add].first;
      added_pairs.push_back(to_add);
      current_value = FastEvalEnergy(sequence);
      std::string current_structure=  BasePairListToStructure1(sequence.length(),current_pairs);
      path.push_back(std::pair<double,std::string>(current_value,current_structure));
      if ( current_value > highest ) {
        highest = current_value;
	highest_structure = current_structure;
      }
    }
  }
  //revert basepair additions and deletions
  //  std::cout<<"PairTable after Partial Path "+PrintPairTable()+"\n";

  for(size_t i=0;i<added_pairs.size();i++){
    pair_table[only_in_base_pairs2[added_pairs[i]].first]=0;
    pair_table[only_in_base_pairs2[added_pairs[i]].second]=0; 
  }
  for(size_t i=0;i<removed_pairs.size();i++){
    pair_table[only_in_base_pairs[removed_pairs[i]].first]=only_in_base_pairs[removed_pairs[i]].second;
    pair_table[only_in_base_pairs[removed_pairs[i]].second]=only_in_base_pairs[removed_pairs[i]].first;
  }

  //  std::cout<<"Restored PairTable "+PrintPairTable()+"\n";
  path.push_back(std::pair<double,std::string>(highest,highest_structure));    
  //return  path;
}

/*
 determine the number of conflicts for each bp in only_in_base_pairs2.Build a group out of the set of bp sharing the lowest # of conflicts 
*/

//std::map<int,std::vector<int> > 
void GetConflictGroup(std::map<int,std::vector<int> > & conflict_group, const std::vector<std::pair<int,int> > & only_in_base_pairs,const std::vector<std::pair<int,int> > & only_in_base_pairs2) {
  conflict_group.clear();
  //bool debug=false;//true;
  std::vector<int> conflicts=std::vector<int>(only_in_base_pairs2.size()); 
  // std::map<int,std::vector<int> >  conflict_group=std::map<int,std::vector<int> >();
  if(only_in_base_pairs2.size()>0){
    size_t minsize=only_in_base_pairs.size();
    for (size_t i=0; i<only_in_base_pairs2.size();i++) {
      conflicts[i]=0;
      for (size_t j=0; j<only_in_base_pairs.size();j++) {
        if(BasePairConflict(only_in_base_pairs2[i],only_in_base_pairs[j])) conflicts[i]++;
      }
      #ifdef _DEBUG_MH_
          std::cout<<"only_in_base_pairs2["+Str((int)i)+"]: ("+Str(only_in_base_pairs2[i].first)+","+Str(only_in_base_pairs2[i].second)+")\n";
          std::cout<<"conflicts[i]:"<<conflicts[i]<<"\n";
      #endif
	  if(conflicts[i]<(int)minsize) minsize=conflicts[i]; 
    }
    #ifdef _DEBUG_MH_
        std::cout<<"minsize:"<<minsize<<"\n";
    #endif
    for (size_t i=0; i<only_in_base_pairs2.size();i++) {
      if(conflicts[i]==(int)minsize){
        conflict_group[i]=std::vector<int>(); 
        for (size_t j=0; j<only_in_base_pairs.size();j++) {
          if(BasePairConflict(only_in_base_pairs2[i],only_in_base_pairs[j])) conflict_group[i].push_back(j);
        }
      }
    }
  }
  else{
    conflict_group[-1]=std::vector<int>();
    for (size_t i=0; i<only_in_base_pairs.size(); i++) conflict_group[-1].push_back(i);
  }
  // return conflict_group;
}



std::vector<std::pair<int,int> > SetDifference(std::vector<std::pair<int,int> > minuend,std::vector<std::pair<int,int> > subtrahend){
  std::vector<std::pair<int,int> > ret=std::vector<std::pair<int,int> > ();
   for (size_t i=0; i<minuend.size(); i++) {
     if (find(subtrahend.begin(),subtrahend.end(),minuend[i]) == subtrahend.end()) ret.push_back(minuend[i]);
   }
   return ret;
}

//std::vector<int> 
void FindBestPartialPathCombination(std::vector<int> & best_combination,int lookahead,std::string sequence,const std::vector<std::pair<int,int> > & backtrack_base,
const std::map<int,std::vector<int> > & conflict_group,const std::vector<std::pair<int,int> > & only_in_base_pairs,
const std::vector<std::pair<int,int> > & only_in_base_pairs2){
 #ifdef _DEBUG_MH_
  std::cout<<"FindBestPartialPathCombination"<<std::endl;
 #endif
  best_combination.clear();
  long int n_combinations = N_take_k( conflict_group.size(),lookahead);
  std::string best_combination_saddle=std::string();
  double best_combination_saddle_energy = INF;//<double>::max();
  //std::cout<<"n_combinations: "+Str((int)n_combinations)<<std::endl;
  best_combination.clear();

 
  std::vector<int> combination= std::vector<int>();  
  std::vector<std::pair<double,std::string> > partial_path=std::vector<std::pair<double,std::string> > ();
  //try all combinations and record combination idx with the best energy.Keep track of saddle point for the entire function all, that is throughout all groups.
  for (long int j=0; j<n_combinations; j++) {
    combination = GetCombination(conflict_group.size(),lookahead,j);
    #ifdef _DEBUG_MH_
    std::cout<<"obtained combination for partial_path"<<std::endl;
    for(size_t r=0;r<combination.size();r++) std::cout<<combination[r];
    std::cout<<std::endl;
    #endif

    DoPartialPath(partial_path,combination,sequence,backtrack_base,conflict_group,only_in_base_pairs,only_in_base_pairs2);
    double current_combination_saddle_energy=partial_path.back().first;
    std::string current_combination_saddle=partial_path.back().second;
     #ifdef _DEBUG_MH_
    std::cout<<"current_combination_saddle "+current_combination_saddle+" "+Str(current_combination_saddle_energy)<<std::endl;
     #endif
    // take the minimum value of the permutations of this group achieved so far.
    if (current_combination_saddle_energy < best_combination_saddle_energy ) {
      best_combination_saddle_energy = current_combination_saddle_energy;
      best_combination_saddle =current_combination_saddle;
      best_combination=combination;
      #ifdef _DEBUG_MH_
        std::cout<<"new best_combination_saddle "+best_combination_saddle+" "+Str(best_combination_saddle_energy)<<std::endl;
        std::cout<<"new best combination ";
        for(size_t r=0;r<best_combination.size();r++) std::cout<<best_combination[r];
        std::cout<<std::endl;
      #endif
    } 
  } 
  // return best_combination;
}

/*
Updates basepairs to reflect the base pairs that have been added to the current structure on the Morgan Higgs path
*/

//std::vector<std::pair<int,int> >  
void UpdateBacktrackBase(std::vector<std::pair<int,int> >& base_pairs,const std::map<int,std::vector<int> > & conflict_group,const std::vector<int> & remove,
const std::vector<std::pair<int,int> > & only_in_base_pairs,const std::vector<std::pair<int,int> > & only_in_base_pairs2){
      for (size_t i=0; i<remove.size(); i++) { 
	int idx = remove[i]-1;
        std::map<int,std::vector<int> >::const_iterator it=conflict_group.begin();
        advance(it,idx);
        std::pair<int,std::vector<int> > conflict=*it;
        int to_add=conflict.first;
        //if(to_add==1) continue;
        std::vector<int> to_remove=conflict.second;  

	for (size_t j=0; j<to_remove.size(); j++) {
	  std::pair<int,int> p=only_in_base_pairs[to_remove[j]];
	  for (size_t k=0; k<base_pairs.size();k++) {
            if(base_pairs[k]==p) {
	      base_pairs.erase(base_pairs.begin()+k);
	      break;
	    }
	  }
	}
	//-1 means that there is nothing to add any more 
	if(to_add!=-1){	  
	  base_pairs.push_back(only_in_base_pairs2[to_add]);
	}
      }
      //  return base_pairs;
}   

/*
bool StackNotMadeFromBases(std::vector<std::pair<int,int> > stack,std::vector<int> bases){
   std::cout<<"StackNotMadeFromBases\n";
   std::cout<<"Base:";
   for(size_t j=0;j<bases.size();j++){
     std::cout<<bases[j]<<" ";
   }
   std::cout<<"\n";
   for(size_t j=0;j<stack.size();j++){
      int base1=stack[j].first;
      int base2=stack[j].second;
      std::cout<<"base1:"+Str(base1)+"\n";
      std::cout<<"base2:"+Str(base2)+"\n";
      //break if base1 was paired
      if(find( bases.begin(), bases.end(),base1)!= bases.end()) return false;
      //break if base2 was paired
      if(find( bases.begin(), bases.end(),base2)!= bases.end()) return false;
    }
   return true;
}
*/

static std::vector<std::pair<double,std::string> > path;
void SeedPath(std::string sequence,std::string src){

path=std::vector<std::pair<double,std::string> >(1,std::pair<double,std::string>(FastEvalEnergy(sequence),src));

}

/**
 *
 */


//std::pair<double,std::string>
std::vector<std::pair<double,std::string> > MorganHiggsStudlaEnergy(std::string sequence,std::string src,std::string tgt,double saddlE,int lookahead,std::string grouping){
  int interrupt=-1;
  bool debug=false;
  #ifdef _DEBUG_MH_
     std::cout<<"MorganHiggsStudlaEnergy with lookahead "+Str(lookahead)+" and grouping "+grouping+"\n";
     std::cout<<"src:\n"+src+":"+Str(FastEvalEnergy(sequence))+"\n"; 
     std::cout<<"tgt:\n"+tgt+":"+Str(EvalEnergy(sequence,tgt))+"\n";
    #endif

  // Obtain the vector of base pairs with values in [1,Node::matrix_size]
  std::vector<std::pair<int,int> > base_pairs  = MakeBasePairList1(src);
  //  std::vector<std::pair<int,int> >
 base_pairs2 = MakeBasePairList1(tgt);
  std::vector<std::pair<int,int> > remove_catalog =SetDifference(base_pairs,base_pairs2);
  std::vector<std::pair<int,int> > add_catalog    =SetDifference(base_pairs2,base_pairs);
  if(MHS_debug) std::cout<<"add_catalog: "+PrintBasePairList(add_catalog)<<endl;
  if(MHS_debug) std::cout<<"remove_catalog: "+PrintBasePairList(remove_catalog)<<endl;
 


  
  //src is the first element in the path
  //  std::vector<std::pair<double,std::string> > path=new std::vector<std::pair<double,std::string> >();


  //path.push_back(std::pair<double,std::string>(FastEvalEnergy(sequence),src));
 
  SeedPath(sequence,src);



  PartialPath partial=PartialPath(base_pairs,base_pairs2,sequence,lookahead,remove_catalog,add_catalog); 
 
  while(partial.add_catalog.size()>0 || partial.remove_catalog.size()>0){    
   
    //treat all combinations in a conflict group
    int prevPathSize=path.size();
    //       MakePairTableFromBasePairs(partial.src,src.size());

       MakePairTable(const_cast<char*>(path.back().second.c_str()));
       if(MHS_debug) cout<<"MorganHiggsEnergy calls FindOptimalPath with pair_table "+PrintPairTable()<<endl;
       partial.FindOptimalPath(path);
       if(MHS_debug) cout<<"pair_table after FindOptimalPath:"+PrintPairTable()<<endl;
       //handle inerrupt
       for(size_t k=prevPathSize;k<path.size();k++){
         if(path[k].first>saddlE && interrupt==-1) interrupt=k;//path.size();   
       }
       //save varialbes between partialpaths
  }
  

  if(interrupt==-1) interrupt=path.size();
  path.push_back(make_pair((double)interrupt,std::string()));
 
  //  cout<<"Obtained path:"<<endl;
  //for(size_t l=0;l<path.size() ;l++){
  //  cout<<path[l].second+" "+Str(path[l].first)<<endl;
  // }
  if(MHS_debug) cout<<PrintPairTable()<<endl;
  return path;
}



//std::pair<double,std::string>
std::vector<std::pair<double,std::string> > MorganHiggsEnergy(std::string sequence,std::string src,std::string tgt,double saddlE,int lookahead,std::string grouping){
  int interrupt=-1;

   #ifdef _DEBUG_MH_
     std::cout<<"MorganHiggsEnergy with lookahead "+Str(lookahead)+" and grouping "+grouping+"\n";
     std::cout<<"src:\n"+src+":"+Str(FastEvalEnergy(sequence))+"\n"; 
     std::cout<<"tgt:\n"+tgt+":"+Str(EvalEnergy(sequence,tgt))+"\n";
      #endif

  // Obtain the vector of base pairs with values in [1,Node::matrix_size]
  std::vector<std::pair<int,int> > base_pairs  = MakeBasePairList1(src);
  //  std::vector<std::pair<int,int> >
 base_pairs2 = MakeBasePairList1(tgt);
  std::vector<std::pair<int,int> > only_in_base_pairs =SetDifference(base_pairs,base_pairs2);
  std::vector<std::pair<int,int> > only_in_base_pairs2 =SetDifference(base_pairs2,base_pairs);

  //Nukleationszentrum: Habe x zusammenhaengende Basenspaare in  only_in_base_pairs2, von denen keine einzige Base in only_in_base_pairs 
  //war. xist die maximale helixgroesse in  base_pairs aber hoechstens 3.
  
 //  std::vector<std::vector<std::pair<int,int> > > front_stacks=ConformationToStacks(base_pairs);
//   std::vector<int> paired_bases_in_front=std::vector<int>(); 
//   int helix_max=0;
//   for(int i=0;i<front_stacks.size();i++){
//     std::vector<std::pair<int,int> > stack= front_stacks[i];
//     for(int j=0;j<stack.size();j++){
//       paired_bases_in_front.push_back(stack[j].first);
//       paired_bases_in_front.push_back(stack[j].second);
//     }
//     if(stack.size()>helix_max) helix_max=stack.size();
//   }
//   if(helix_max>3) helix_max=3;
//   std::cout<<"helix_max:"+Str(helix_max )+"\n";


//   std::cout<<"# elts in only_in_base_pairs2:"+Str((int)only_in_base_pairs2.size())+" \n";
//   std::vector<std::vector<std::pair<int,int> > > new_stacks=ConformationToStacks(only_in_base_pairs2);
//   std::cout<<"# new stacks:"+Str((int)new_stacks.size())+" \n";
//   bool have_nucleation=false;
//   for(int i=0;i<new_stacks.size();i++){
//     std::vector<std::pair<int,int> > stack= new_stacks[i];
//     if(stack.size()<helix_max) continue;
//     //the stack is large enough. Ensure that none of its bases are in a base pair in base_pairs.
//     if(StackNotMadeFromBases(stack,paired_bases_in_front)){
//       std::cout<<"Have nucleation\n";
//       have_nucleation=true;
//       break;
//     }
//   }
//   if(!have_nucleation){
//     std::cout<<"Cannot fold form src to tgt because there is no nucleation center\n";
//     std::cout<<src+"\n";
//     std::cout<<tgt+"\n";
//     std::vector<std::pair<double,std::string> > path= std::vector<std::pair<double,std::string> >();
//     path.push_back(std::pair<double,std::string>(1000,src));
//     return path;
//   }
  
  //src is the first element in the path
  std::vector<std::pair<double,std::string> > path=std::vector<std::pair<double,std::string> >();
  path.push_back(std::pair<double,std::string>(FastEvalEnergy(sequence),src));
  std::vector<std::pair<double,std::string> > partial_path= std::vector<std::pair<double,std::string> >();
  std::map<int,std::vector<int> > conflict_group= std::map<int,std::vector<int> >();
  std::vector<int> combination= std::vector<int>();
  while(only_in_base_pairs2.size()>0 || only_in_base_pairs.size()>0){
     
    #ifdef _DEBUG_MH_
       std::cout<<"base_pairs:\n"<<BasePairListToStructure1(sequence.length(),base_pairs)<<"\n";
       std::cout<<"base_pairs2:\n"<<BasePairListToStructure1(sequence.length(),base_pairs2)<<"\n";
       std::cout<<"only_in_base_pairs:\n"<<BasePairListToStructure1(sequence.length(),only_in_base_pairs)<<"\n";
       std::cout<<"only_in_base_pairs2:\n"<<BasePairListToStructure1(sequence.length(),only_in_base_pairs2)<<"\n";
    #endif
      
       GetConflictGroup(conflict_group,only_in_base_pairs,only_in_base_pairs2); 
    #ifdef _DEBUG_MH_
       std::cout<<"got conflict_group"<<"\n";
       std::cout<<"PrintConflictGroup\n";
       std::cout<<PrintConflictGroup(conflict_group,only_in_base_pairs,only_in_base_pairs2)+"\n";
    #endif
    //treat all combinations in a conflict group
    while(conflict_group.size()>0){
      #ifdef _DEBUG_MH_
          std::cout<<"conflict_group.size "+Str((int)conflict_group.size())<<std::endl;
          std::cout<<"Look for a path combination"<<std::endl;
      #endif
	
      FindBestPartialPathCombination(combination,lookahead,sequence,base_pairs,conflict_group,only_in_base_pairs,only_in_base_pairs2);  
      #ifdef _DEBUG_MH_
          std::cout<<"PrintCombination()"<<std::endl;
          std::cout<<PrintCombination(combination)<<std::endl;
          std::cout<<"take its partial path"<<std::endl;
      #endif 

       DoPartialPath(partial_path,combination,sequence,base_pairs,conflict_group,only_in_base_pairs,only_in_base_pairs2); 
      
      //Stop one element short as to not push group saddle back here. It is done right before MorganHiggsEnergy returns!
      #ifdef _DEBUG_MH_
          std::cout<<"Print partial path"<<std::endl;
      #endif
      for(size_t k=0;k<partial_path.size()-1;k++) {
	 #ifdef _DEBUG_MH_
	std::cout<<partial_path[k].second+":"+Str(partial_path[k].first)<<std::endl;
         #endif
         path.push_back(partial_path[k]);
	 //XXXXXXXXXXXXXX
	 //can only interrupt once
         if(partial_path[k].first>saddlE && interrupt==-1) interrupt=path.size();   
      }


      //update the elements in base_pairs.
      #ifdef _DEBUG_MH_
          std::cout<<"UpdateBacktrackBase"<<std::endl;
          std::cout<<"old backtrackbase:\n"<<BasePairListToStructure1(sequence.length(),base_pairs)<<std::endl;
      #endif
      UpdateBacktrackBase(base_pairs,conflict_group,combination,only_in_base_pairs,only_in_base_pairs2);
      MakePairTableFromBasePairs(base_pairs,src.size());
      #ifdef _DEBUG_MH_
          std::cout<<"new backtrackbase:\n"<<BasePairListToStructure1(sequence.length(),base_pairs)<<std::endl;
          //find new conflict group after each constructed partial path
          std::cout<<"grouping is _"+grouping+"_"<<std::endl; 
      #endif
      if(grouping=="regroup") break;
      else if(grouping=="standard"){

          #ifdef _DEBUG_MH_REGROUP_
             std::cout<<"Regroup ConflictGroup"<<std::endl;
             std::cout<<PrintConflictGroup(conflict_group,only_in_base_pairs,only_in_base_pairs2)+"\n";
             //std::cout<<PrintConflictGroup(conflict_group)<<std::endl;
	  #endif 
          sort(combination.begin(),combination.end());
          for (size_t k=0; k<combination.size(); k++) { 
          //delete elements in the combination, i.e. those that have not been used to construct the partial path
	    //add 1 to k, as combinations start at 1, not at 0
            //remove the ones that _were_ in the combination
	    //int k=it->first+1;
            // combination 1 2
            //conflict group
            //0: 0 1
            //1: 0 1
            //3: 1 2
            //need:
	    //1 -> 0: 0 1
	    //2 -> 1: 0 1
            std::map<int,std::vector<int> >::iterator it=conflict_group.begin();
            //the -k compensates for the deletions from conflict_group in previous iterations
            advance(it,combination[k]-1-k);       
	    #ifdef _DEBUG_MH_REGROUP_
                //the index of the base pair in only_in_base_pairs2==oibp2 that k refers to in the combination and that hence has been added to the path
                size_t idx_oibp2=it->first;
                std::cout<<"before erasing "+Str(idx_oibp2)<<std::endl;
                std::cout<<PrintConflictGroup(conflict_group,only_in_base_pairs,only_in_base_pairs2)+"\n";
                //std::cout<<PrintConflictGroup(conflict_group)<<std::endl;
	     #endif 
             conflict_group.erase(it);
             #ifdef _DEBUG_MH_REGROUP_
                std::cout<<"after erasing "+Str(idx_oibp2)<<std::endl;
                std::cout<<PrintConflictGroup(conflict_group,only_in_base_pairs,only_in_base_pairs2)+"\n";
                //std::cout<<PrintConflictGroup(conflict_group)<<std::endl;
             #endif 
	  }
        #ifdef _DEBUG_MH_
            std::cout<<"New ConflictGroup"<<std::endl;
            std::cout<<PrintConflictGroup(conflict_group,only_in_base_pairs,only_in_base_pairs2)+"\n";
            //std::cout<<PrintConflictGroup(conflict_group)<<std::endl;
        #endif 
      }
    }
    #ifdef _DEBUG_MH_
        std::cout<<"recalculate only_in_base_pairs"<<std::endl;
    #endif 
    only_in_base_pairs=SetDifference(base_pairs,base_pairs2);
    #ifdef _DEBUG_MH_
        std::cout<<"recalculate only_in_base_pairs2"<<std::endl;
    #endif 
    only_in_base_pairs2=SetDifference(base_pairs2,base_pairs);
  }
  
//   //the global saddle, the max across the final path
//   std::pair<double,std::string> saddle=path[0];
//   //std::cout<<"MH path"<<std::endl;
//   for(size_t i=1;i<path.size();i++){
//     if(path[i].first>saddle.first) saddle=path[i];
//     //std::cout<<path[i].second+" "+Str(path[i].first)<<std::endl;
//   }
//   path.push_back(saddle);
  
  if(interrupt==-1) interrupt=path.size();
  path.push_back(make_pair((double)interrupt,std::string()));
  return path;
}


std::string PrintConflictGroup(std::vector<std::vector<int> > conflict_group,std::vector<std::pair<int,int> > only_in_base_pairs,
std::vector<std::pair<int,int> > only_in_base_pairs2){
  std::string s;
  for(std::vector<std::vector<int> >::iterator it=conflict_group.begin();it!=conflict_group.end();it++){
        std::vector<int> v=*it;
        s+=Str(v[0])+" "+PrintBasePair(only_in_base_pairs2[v[0]])+":";
        for(size_t j=1;j<v.size();j++) s+=Str(v[0])+" "+PrintBasePair(only_in_base_pairs[v[j]])+" ";
        s+="\n";
      }
  return s;

  //  int to_add=conflict.first;
    
  //  std::vector<int> to_remove=conflict.second;
}
//ORIGINAL
std::string PrintConflictGroup(std::map<int,std::vector<int> > conflict_group,std::vector<std::pair<int,int> > only_in_base_pairs,
std::vector<std::pair<int,int> > only_in_base_pairs2){
  std::string s;
  for(std::map<int,std::vector<int> >::iterator it=conflict_group.begin();it!=conflict_group.end();it++){
        s+=PrintBasePair(only_in_base_pairs2[it->first])+":";
        std::vector<int> v=it->second;
        for(size_t j=0;j<v.size();j++) s+=PrintBasePair(only_in_base_pairs[v[j]])+" ";
        s+="\n";
      }
  return s;

  //  int to_add=conflict.first;
    
  //  std::vector<int> to_remove=conflict.second;
}

std::string PrintCombination(std::vector<int> combination){
  std::string s;
  for(std::vector<int>::iterator it=combination.begin();it!=combination.end();it++){
        s+=Str(*it)+" ";
      }
  s+="\n";
  return s;

}

void MakePairTableFromBasePairs(const std::vector<std::pair<int,int> >& base_pairs,int length){
 for(int i=1;i<=length;i++){
   pair_table[i]=0;
 }
 for(size_t i=0;i<base_pairs.size();i++){  
   pair_table[base_pairs[i].first]=base_pairs[i].second;
   pair_table[base_pairs[i].second]=base_pairs[i].first;
  }

}


#ifdef _TEST_MORGANHIGGS_
extern "C" {
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"
#include "pair_mat.h"
}
//extern 
short *S;
//extern 
short *S1;
//extern 
short *pair_table;
// End of file
int main(int argc, char *argv[]) {

  double T=-2.2;
  double mod= fmod(T, 1);
  std::cout<<mod<<std::endl;

  OptionS* OptS;
  OptS = decodeCML(argc, argv);
  double saddlE=10000;
  std::cout<<"Enter sequence, start structure and target structure, separated by new lines"<<std::endl;
  std::string sequence;
  std::string src;
  std::string tgt;

  std::cin>>sequence;
  std::cin>>src;
  std::cin>>tgt;
  OptS->transcribed=src.size();
  InitializeEnergyModel(OptS,sequence);
  //InitializeViennaRNA(sequence,OptS->dangle,src.size());
  MakePairTable(const_cast<char*>(src.c_str()));
  std::cout<<"new:"+PrintPairTable()<<std::endl;

  std::cout<<"S "<<S[0]<<std::endl;
  std::cout<<"S1 "<<S[1]<<std::endl;

  std::vector<std::pair<double,std::string> > path;
  if(OptS->barrier_heuristic=='M'){
      path=MorganHiggsEnergy(sequence,src,tgt,saddlE,OptS->lookahead,OptS->grouping);
  }
  else if(OptS->barrier_heuristic=='S'){
     path=MorganHiggsStudlaEnergy(sequence,src,tgt,saddlE,OptS->lookahead,OptS->grouping);
  }
  for(size_t l=0;l<path.size() ;l++){
    cout<<path[l].second+" "+Str(path[l].first)<<endl;
  }
  return 0;
}
#endif
