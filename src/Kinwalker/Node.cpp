/*
  Last changed Time-stamp: <2007-07-10 19:27:22 xtof>
  $Id: Node.cpp,v 1.35 2007/11/03 16:45:58 Kinwalker Exp $
*/

#include <algorithm>
#include <iostream>
#include <iterator>
#include <utility>
#include "Node.h"

#include "template_utils.c"

extern "C" {
#include "findpath.h"
}

#define MS_PER_TIME_UNIT .0001
#define TIME_VS_DELTAG_DY_DX (8.0/11.0)
#define EPSILON .00000000001

#define MIN_ENERGY_DIFF .01

double BARRIER_TOO_HIGH=10000;//std::numeric_limits<int>::max();
// class variables
int Node::verbose;
int Node::lookahead;
int Node::matrix_size;
std::string Node::grouping;
std::string Node::sequence;
std::string Node::front_structure;
double Node::front_energy;
std::string Node::mfe_structure;
double Node::mfe;
double Node::energy_barrier;
double Node::max_barrier;
std::vector<std::vector<bool> > Node::front =std::vector<std::vector<bool> >();
std::vector<Node*> Node::front_extrema = std::vector<Node*>() ;
std::vector<Node*> Node::extrema = std::vector<Node*>() ;
std::vector<double> Node::extension_cost = std::vector<double>() ;


std::string Node::constraint_string;
int Node::ineligible_count=0;
int Node::transcribed;
int Node::initially_transcribed;
std::string Node::front_trajectory_ps;
bool Node::print_front_trajectory=false;
int Node::interrupt_trajectory;
int Node::windowsize;

bool Node::transcribed_front_maximal=false;
double Node::transcription_rate;
double Node::time;
OptionS* Node::OptS=NULL;
std::vector<std::string>  Node::trajectory= std::vector<std::string>();
std::vector<std::pair<double,std::string> > Node::path=std::vector<std::pair<double,std::string> > ();

short *S;
short *S1;
short *pair_table;



struct BasePairConflict :public std::binary_function<std::pair<int,int>,std::pair<int,int>,bool>{
   bool operator()(std::pair<int,int> p1,std::pair<int,int> p2) const{
     if(p1.first==p2.first && p1.second==p2.second) return false;
     if(p1.first==p2.first || p1.first==p2.second || p1.second==p2.first || p1.second==p2.second) return true;
     if(p1.first<=p2.first && p1.second>= p2.first && p1.second<=p2.second) return true;
     if(p1.first>=p2.first && p1.first<=p2.second && p1.second>=p2.second)  return true;    
     else return false;
   }
};


double Node::DeltaGToTimeUnits(double delta_G){
  //Cout("delta_G:"+Str(delta_G));
  double units= pow(10.0,TIME_VS_DELTAG_DY_DX*delta_G);

  //Cout("units:"+Str(units));
  return units;
}

double Node::TimeUnitsToDeltaG(double units){
  return log10(units)/TIME_VS_DELTAG_DY_DX;
}

double Node::TimeToDeltaG(double time){
  if(time<=EPSILON) return 0.0;
  double units=time*1000/MS_PER_TIME_UNIT;
  double delta_G= Node::TimeUnitsToDeltaG(units);
  if(delta_G<0){
    Cout("Error in TimeToDelta G. Delta G is "+Str(delta_G)+" for time "+Str(time));
    exit(0);
  } 
  return delta_G;
}

double Node::TimePassedFolding(double delta_G){
  double units=Node::DeltaGToTimeUnits(delta_G);
  //  Cout("units:"+Str(units)+"\n");
  double msec=units*MS_PER_TIME_UNIT;
  // Cout("msec:"+Str(msec)+"\n");
  if(verbose>=2) Cout(Str(msec*.001)+" secs passed folding from dG "+Str(delta_G)+"\n"); 
  return msec*.001;
}



 void  Node::ProcessOptions(OptionS* OptS){ 
  // Allocate memory
  //obtain input sequence
  Node::OptS=OptS;
  if (OptS->testseq) Node::sequence = "ACAGGUUCGCCUGUGUUGCGAACCUGCGGGUUCG";
  else std::cin >> Node::sequence;


  Node::print_front_trajectory=OptS->printfront;
  Node::grouping = std::string(OptS->grouping); 
  Node::verbose = OptS->verbose;
  Node::lookahead = OptS->lookahead;
  Node::matrix_size = Node::sequence.size();
  if(OptS->windowsize<=0) Node::windowsize= Node::sequence.size();
  else {
    Node::windowsize=OptS->windowsize;
    if(Node::windowsize>(int)Node::sequence.size()) Node::windowsize=Node::sequence.size();
  }
  Node::transcribed=OptS->transcribed;
  Node::initially_transcribed=Node::transcribed;
  Node::interrupt_trajectory=OptS->interrupt_trajectory;
  if(Node::transcribed<0) Node::transcribed=Node::matrix_size;

  Node::front_structure=std::string(Node::transcribed,'.');

  if(OptS->init_structure) {
    Cout("Please enter starting structure\n");
    std::cin >> Node::front_structure;
  }
  // fold_constrained=OptS->fold_constrained;
//   if(fold_constrained) {
//     Cout("Please enter folding constraints\n");
//     std::cin >> Node::constraint_string;
//   }

  Node::max_barrier=0.0; 
  Node::energy_barrier=0.0; 
  if(Node::verbose>=1) std::cout<<"Initial number of bases transcribed:"+Str(Node::transcribed)+"\n";
  Node::transcription_rate=(double)OptS->transcription_rate;
  if(Node::verbose>=1) std::cout<<"Transcription rate:"+Str(Node::transcription_rate)+"\n";
  Node::time=0.0;
  Node::CalculateMaxBarrier();
 
 
}

void Node::CalculateMaxBarrier(){
  double t_transcribe=(double)1.0/(double)Node::transcription_rate;
  Node::max_barrier=Node::TimeToDeltaG(t_transcribe);
}




/**
 * Constructor
 */
Node::Node(int i, int j, double energy)
{
  this->i = i;
  this->j = j;
  this->free_energy = energy;
  this->is_eligible=true;
  this->is_included=false;
}


std::string Node::FoldingPathToString(std::vector<std::pair<double,std::string> > path){
  std::string ret="";
  for(size_t i=0;i<path.size();i++)
    ret+=path[i].second+" "+Str(path[i].first)+" "+Str(Node::energy_barrier)+"\n";
  return ret;
}

void Node::GetMorganHiggsPath(std::string target){

  int t=Node::transcribed;

  Node::path = MorganHiggsEnergy(Node::sequence.substr(0,t),
                                 Node::front_structure.substr(0,t),
                                 target.substr(0,t),
                                 Node::front_energy+Node::energy_barrier,
                                 Node::lookahead,Node::grouping);//,Node::interrupt_trajectory);
}

void Node::GetMorganHiggsStudlaPath(std::string target){

  int t=Node::transcribed;

  Node::path = MorganHiggsStudlaEnergy(Node::sequence.substr(0,t),
                                 Node::front_structure.substr(0,t),
                                 target.substr(0,t),
                                 Node::front_energy+Node::energy_barrier,
                                       Node::lookahead,Node::grouping);//,Node::interrupt_trajectory);
}

/**
counts from 0.
*/
int Node::FindPathMinimum(int last_idx_within_barrier){

   double min=INF;
   int min_idx=-1;
   for(size_t i=0;i<=(size_t)last_idx_within_barrier;i++){
     if(Node::path[i].first<min){
       min= Node::path[i].first;
       min_idx=i;
     }
   }
   //Cout("FindPathMinimum "+Str(Node::path[min_idx].first)+"\n");
   return min_idx;
}


/**
counts from 0.
*/
int Node::FindPathSaddle(int min_idx){
  bool debug=false;
  if(debug) Cout("FindPathSaddle got min_idx"+Str(min_idx)+"\n");
  double max=-INF;
   int max_idx=-1;
   for(int i=0;i<=min_idx;i++){
      if(debug) Cout("i "+Str(i)+" "+Str(Node::path[i].first)+" max "+Str(max)+"\n");
     if(Node::path[i].first>max){
       
       max=Node:: path[i].first;
       max_idx=i;
       if(debug)Cout("new max "+Str(max)+" max_idx "+Str(max_idx)+"\n");
     }
   }
   if(debug) Cout("FindPathSaddle idx"+Str(Node:: path[max_idx].first)+"\n");
   return max_idx;
}


void Node::GetSaddleFromPath(std::pair<double,std::string> & saddle,std::pair<double,std::string> & final_structure){  
  //the last element is the dummy element, subtract it here,subtract another 1 for counting from 0
  int path_size=(int)path.size()-2;
  if(Node::interrupt_trajectory){ 
    int last_idx_of_partial_path=(int)path.back().first;
    int min_idx;
    while(true){    
      //Cout("last_idx_of_partial_path "+Str(last_idx_of_partial_path)+" path.size "+Str(path_size));
      min_idx=FindPathMinimum(last_idx_of_partial_path);
      //if this energy of the min within reach is not low enough, we take the next min
      if(path[min_idx].first+MIN_ENERGY_DIFF<Node::front_energy) break;
      while(last_idx_of_partial_path<=path_size && path[last_idx_of_partial_path].first+MIN_ENERGY_DIFF>=Node::front_energy) last_idx_of_partial_path++;
      //if we have increased  to the size of the path, it means there is no structure that is lower, in that case retrun the front structure,
      //as this trajectory is worthless, even when the energy barrier is higher
      if(last_idx_of_partial_path>path_size){
        final_structure=std::make_pair(Node::front_energy,Node::front_structure);
        return;
      }
    }
    int saddle_idx = FindPathSaddle(min_idx);
    saddle=path[saddle_idx];
    final_structure=path[min_idx];
  }
  else{
      int saddle_idx=FindPathSaddle(path_size);
      saddle=path[saddle_idx];
      final_structure=path[path_size];
  }
  //  if(verbose>=4) Cout("GetSaddleFromPath obtains saddle and final structue\n"+saddle.second+"\n"+final_structure.second);
}

/**
Returns a pair consisting of the barrier height to overcome to fold into the new structure from the front as first element
and the new structure as the second element.
*/



void Node::CalculateFoldingPath(Node* extremum,std::string integrated_structure){
  if(verbose>=3) std::cout<<"Folding "<<OptS->barrier_heuristic<<std::endl;
    MakePairTableFromFrontStructure();
    if(OptS->barrier_heuristic=='M')
      GetMorganHiggsPath(integrated_structure);
    else if(OptS->barrier_heuristic=='S') 
      GetMorganHiggsStudlaPath(integrated_structure);
    else { /* other heuristic */
      path_t *p;
      int p_len = 0;
      double maxE = -INF;//std::numeric_limits<double>::max();
      int maxE_idx = 0;
      int t = Node::transcribed;
      std::vector<std::pair<double,std::string> > v;
      p = get_path(const_cast<char*>(sequence.substr(0,t).c_str()),
                   const_cast<char*>(Node::front_structure.c_str()),
                   const_cast<char*>(integrated_structure.substr(0,t).c_str()),
                   Node::OptS->maxkeep);
      for (int i=0; p[i].s != NULL; i++) {
        // memorize idx of structure with highest energy seen so far
        if(p[i].en > maxE){
          maxE = p[i].en;
          maxE_idx = i+1;
        }
        v.push_back(std::make_pair(p[i].en, p[i].s));
      }
      // add dummy entry with idx of structure with highest energy
      v.push_back(std::make_pair(maxE_idx, ""));       
      /* clean up space of path */
      for (int i=0; p[i].s != NULL; i++) free(p[i].s);
      free(p);
      Node::path = v;

    }

    if(verbose>=3) {
      Cout(Node::front_structure+"\n");
      Cout(integrated_structure+"\n");
      Cout(Node::FoldingPathToString(Node::path));
    }
}

void Node::CalculatePathAndSaddle(Node* extremum,std::string integrated_structure,std::pair<double,std::string> & saddle,std::pair<double,std::string> & final_structure){

  //get path and saddle for each heuristic separately. Keep the one that produced the lowest saddle
  char heuristic=OptS->barrier_heuristic;

  OptS->barrier_heuristic='M';
  Node::CalculateFoldingPath(extremum,integrated_structure);
  Node::GetSaddleFromPath(saddle,final_structure);
  std::pair<double,std::string> new_saddle=std::pair<double,std::string> ();
  std::pair<double,std::string> new_final_structure=std::pair<double,std::string> ();
  OptS->barrier_heuristic='S';

  Node::CalculateFoldingPath(extremum,integrated_structure);
  Node::GetSaddleFromPath(new_saddle,new_final_structure);
  if(new_saddle.first<saddle.first && StructureIsFrontExtension(new_final_structure)){
    saddle=new_saddle;
    final_structure=new_final_structure;
  }
  OptS->barrier_heuristic='F';
  Node::CalculateFoldingPath(extremum,integrated_structure);
  Node::GetSaddleFromPath(new_saddle,new_final_structure);
  if(new_saddle.first<saddle.first && StructureIsFrontExtension(new_final_structure)){
    saddle=new_saddle;
    final_structure=new_final_structure;
  }

  OptS->barrier_heuristic=heuristic;
}


std::string Node::CombineFrontAndNode(std::string node_substructure){
  bool debug=false;
  if(debug) Cout("CombineFrontAndNode\n front "+Node::front_structure+"\n node  "+node_substructure+"\n");
  
  std::vector<std::pair<int,int> >bp_front=MakeBasePairList1(Node::front_structure);
  std::vector<std::pair<int,int> >bp_node=MakeBasePairList1(node_substructure);
  for(size_t i=0;i<bp_front.size();i++){
    std::pair<int,int> bp=bp_front[i];
    //do nothing if the base pair is already part of basepairs
    if(find(bp_node.begin(),bp_node.end(),bp)!=bp_node.end()) continue;
    if(!Conflict(bp_node,bp)) {
        if(debug) Cout("added bp "+PrintBasePair(bp)+"\n");
        node_substructure[bp.first-1]='(';
        node_substructure[bp.second-1]=')';
      }
    }
  if(debug) Cout("combined "+node_substructure);
    return node_substructure;
}


bool Node::IsWithinWindow(){
  if(this->j-this->i<=Node::windowsize) return true;
  else return false;
}

/**
 * Add the next extremum that does not conflict with the front and is
 * within the limits imposed by energy_barrier.
 */
bool
Node::FindExtremum(){
  std::pair<double,std::string> saddle= std::pair<double,std::string>();
  std::pair<double,std::string> final_structure= std::pair<double,std::string>();
  Node::transcribed_front_maximal=false;
  //return if no extrema are transcribed
  if(find_if(extrema.begin(),extrema.end(),std::mem_fun(&Node::IsTranscribed))==extrema.end()){
    Node::transcribed_front_maximal=true;
    return false;
  }
  extension_cost.clear();
  MakePairTableFromFrontStructure();
  std::vector<std::pair<int,int> >bp_front=MakeBasePairList1(front_structure);
  for (size_t i=0; i<extrema.size(); i++) {
    //make sure the value from the last iteration is not reused
    std::string combined_structure="";
    if(!extrema[i]->IsTranscribed() || !extrema[i]->IsEligible() || (!AllTranscribed() && !extrema[i]->IsWithinWindow()) ) {
      //if(verbose>=2) Cout("Not transcribed.\n");
      extension_cost.push_back(BARRIER_TOO_HIGH);
      continue;
    }
    if(verbose>=4) Cout("Try extending to Node "+Print(extrema[i])+".\n");
    if(verbose>=3) Cout("Front is\n"+Node::front_structure+":"+Str(Node::front_energy)+"\n");
    //try adding terminal basepairs if backtracked structure doesn't add anything to front
    bool debug=false;
    double barrier=BARRIER_TOO_HIGH;
    std::string node_substructure=BacktrackNode(extrema[i]);
    if(debug)    Cout("node_substructure\n"+node_substructure);
    std::string integrated_structure=CombineFrontAndNode(node_substructure);
    bool skip=false;
    if(!Node::interrupt_trajectory){
      if(integrated_structure==Node::front_structure || Evaluate(integrated_structure)>Node::front_energy) skip=true;
    }

    if(verbose>=3) Cout("Combination of Node and Front\n"+integrated_structure+"\n");
    if(!skip && Node::front_structure!=integrated_structure){
      if(OptS->barrier_heuristic=='A') CalculatePathAndSaddle(extrema[i],integrated_structure,saddle,final_structure);
      else{
         Node::CalculateFoldingPath(extrema[i],integrated_structure);
         Node::GetSaddleFromPath(saddle,final_structure);
      }
      //we can only accept the final structure if it is different from the front structure and has a better energy
      if(StructureIsFrontExtension(final_structure)) {
        barrier=saddle.first-Node::front_energy;
        combined_structure=final_structure.second;
      }
    }
     if(debug) Cout("Current combined structure and barrier\n"+combined_structure+" "+Str(barrier));
    //needs to be set to a value that gives this one a second chance below once all is transcribed
    //else combined_structure=Node::front_structure;
    //only if everything is transcribed and we are moving towards the mfE
    //done by the function CalculateFoldingPath
    skip=false;

    if(!Node::interrupt_trajectory){
      if(node_substructure==Node::front_structure || Evaluate(node_substructure)>Node::front_energy) skip=true;
    }

    if(!skip && barrier>Node::energy_barrier  && combined_structure!=node_substructure){
      if(OptS->barrier_heuristic=='A') CalculatePathAndSaddle(extrema[i],node_substructure,saddle,final_structure);
      else{
        Node::CalculateFoldingPath(extrema[i],node_substructure);
        Node::GetSaddleFromPath(saddle,final_structure);
      }
      if(StructureIsFrontExtension(final_structure)) {
        double barrier2=saddle.first-Node::front_energy;
        if(barrier2<barrier){
          barrier=barrier2;
          combined_structure=final_structure.second;
          if(debug) Cout("Have valid front extension with\n"+combined_structure+" "+Str(barrier));
        }
      }
    }
    if(debug) Cout("Current combined structure and barrier\n"+combined_structure+" "+Str(barrier));
    if(verbose>=3) Cout("Actually obtained combination:\n"+combined_structure+"\n");
    extension_cost.push_back(barrier);
    

    if(verbose>=2) {
      std::cout<<"Node "+Print(extrema[i])+" has extension cost "+Str(extension_cost.back())+" ("+Str(Node::energy_barrier)+" "+Str(Node::transcribed)+")\n";
    }
    if ( barrier <= Node::energy_barrier ) {
       //1.First structure if that is suitable.
       //2.Otherwise seconds structure, if suitable, or we won't get here
       extrema[i]->AddToFront(combined_structure,barrier);
       if(find_if(extrema.begin()+i,extrema.end(),std::mem_fun(&Node::IsTranscribed))==extrema.end()) Node::transcribed_front_maximal=true;
       return true;
     }
     //Everything is eligible 
     else if(barrier > Node::max_barrier) extrema[i]->SetIneligible();
     else{
       if(verbose>=3)std::cout<<"Extremum"+Print(extrema[i])+" has extension cost "+Str(barrier)+" while "+Str(Node::energy_barrier)+" is allowed and max is "+Str(Node::max_barrier)+"\n";
     }
  }
  return (false);
}

bool Node::StructureIsFrontExtension( const std::pair<double,std::string> & struc){
  if(struc.second!=Node::front_structure && struc.first<Node::front_energy) return true;
  else return false;
}


static int front_count=1;
void
Node::AddToFront(std::string new_front_structure,double barrier){
   this->Reeligify(Node::front_structure,new_front_structure);
    double new_front_energy =Evaluate(new_front_structure);
    if ( verbose >= 1 ) {
         std::cout<<"Add Node ("+Print(this)+" "+Str(extension_cost.back())+" ("+Str(Node::energy_barrier)+" "+Str(Node::transcribed)+")\n";        
         std::cout<<"src:"+Node::front_structure+" "+Str(Node::front_energy)+"\n";
         std::cout<<"tgt:"+new_front_structure+" "+Str(new_front_energy)+"\n";
        }
        
        //add the extremum to the front
        for (int i=this->i; i<=this->j; i++) {
          for (int j=i; j<=this->j; j++) front[i][j] = true;
        }
        front_extrema.push_back(this);
        Node::front_structure=new_front_structure;
        Node::front_energy=new_front_energy;

        if(barrier>0.0) Node::IncreaseTime(TimePassedFolding(barrier));
        std::string trajectory_entry=std::string();
        trajectory_entry+=Node::front_structure+" ";
        trajectory_entry+=Str(Node::front_energy)+" ";
        trajectory_entry+=Str(Node::time)+" ";
        trajectory_entry+=Str(barrier)+" ";
        trajectory_entry+=Str(Node::energy_barrier)+" ";
        trajectory_entry+=Str(Node::transcribed)+" ";
        trajectory.push_back(trajectory_entry);
   
        //and remove it and the extrema it dominates from the extrema vector 
        int count=0;
        for(size_t i=0;i<extrema.size();i++){
          if(!extrema[i]->is_included && extrema[i]->i>=this->i && extrema[i]->j<=this->j){
            if(extrema[i]->IsMfE()) continue;
            extrema[i]->SetIncluded(true);
            extrema[i]->SetIneligible();
            count++;
          }
        }
        
        //print new front
        if ( verbose >= 1 ) {
          std::cout<<"Path to new front:\n";
          std::cout<<FoldingPathToString(this->path)<<std::endl;
          Cout("New Front\n");
          std::cout<<front_structure+" "+Str(front_energy)<<std::endl; 
        }

        if ( verbose >= 5 ) Cout(Node::PrintFront()+"\n");
        //make all extrema eligible again.
        if(Node::print_front_trajectory) {
          std::vector<std::pair<int,int> > ext= std::vector<std::pair<int,int> >();
          for(size_t k=0;k<front_extrema.size();k++)ext.push_back(std::pair<int,int> (front_extrema[k]->i,front_extrema[k]->j));
          Node::front_trajectory_ps+=PSFrontPlot(Node::sequence,ext);

          //          /*
         Node::front_trajectory_ps=PSFrontPlot(Node::sequence,ext);
         Node::front_trajectory_ps+="end \n";
         std::string filename="front_trajectory"+Str(front_count)+".ps";
         std::ofstream outfile(filename.c_str());
         outfile << Node::front_trajectory_ps;
         outfile.close();
         front_count++;
         //*/

        }
}


/**
 * Extends the front by as many non-conflicting extrema as possible
 * within the limits imposed by Node::energy_barrier.
 */
//void
bool
Node::ExtendFront()
{
  bool extended= FindExtremum();
  return extended;
}


/**
 * Find the local extrema in the max_matching matrix and sort them
 * comparing with Node::LessThan
*/

extern "C" {
void export_fold_arrays(int **f5_p, int **c_p, int **fML_p, int **fM1_p,
                        int **indx_p, char **ptype_p);
}
void
Node::FindLocalExtrema()
{
  int *f5, *c, *fML, *fM1, *indx; char *ptype;
  export_fold_arrays(&f5,&c,&fML,&fM1,&indx,&ptype);
  if(verbose>=2) Cout("#extrema: "+Str((int)extrema.size())+"\n");
  int n = matrix_size;
  for (int i=1; i<=n-TURN-1; i++) {
    for (int j=i+TURN+1; j<=n; j++) {
      double val = c[indx[j]+i];
      //i and j pair AND no lonely pair
      if(val < INF) {
        //Node node=Node(i, j,val);
        extrema.push_back(new Node(i, j,val));
        if(verbose>=2) Cout(Print(extrema.back())+" "+Str(val)+"\n");
      }
    }
  }

  double energy = f5[Node::sequence.size()];
  // (1,n) is in general not a pair. You can't include it in the list
  //Node node2=Node(1,n,energy);
  extrema.push_back(new Node(1,n,energy));//& node2);
  if(verbose>=3) Cout("#extrema: "+Str((int)Node::extrema.size())+"\n");
  stable_sort(extrema.begin(), extrema.end(), Node::LessThan);  
}


/**
 * Counts the number of squares by which the Node increases the
 * current front.
 */
int
Node::FrontIncrease()
{
  //WE COUNT substructures FROM  1, but front starts at front[0][0]
  //  std::cout<<"FrontIncrease\n";
  int count = 0;
  for (int i=this->i; i<=this->j; i++) {
    for (int j=i; j<=this->j; j++) {
      // std::cout<<"i: "+Str(i)+" j:"+Str(j)+"\n";
      if ( !front[i-1][j-1] ) count++;
    }
  }
  return count;
}


bool Node::TrulyContainsFront(){
 for (int i=1; i<=Node::matrix_size; i++) {
      if (front[i][this->j] ) return false;
 }
 for (int j=1; j<=Node::matrix_size; j++) {
      if (front[this->i][j] ) return false;
 }
 return true;
}

/**
 * Determines whether Node::front increases when this is added to
 * it. Used to filter extrema, that will not extend the front.
 */
bool
Node::FrontIncreases()
{
  return ( FrontIncrease() > 0 );
}

/**
 * Comparison for STL's sort. Comparison criteria: Distance from main
 * diagonal, distance to left/right edge of the diagonal a Node is on,
 * 5' before 3'.
 */
inline bool
Node::LessThan(Node* n1,Node* n2)
{
  // sort left/right edge of diagonal first
  if ( n1->j-n1->i == n2->j-n2->i ) {
    // Distance to left/right edge of diagonal
    int edge_dist  = std::min(n1->i-1,Node::matrix_size+1-n1->j);
    int edge_dist2 = std::min(n2->i-1,Node::matrix_size+1-n2->j);

    if ( edge_dist == edge_dist2 ) {
      // Take the point on the 5' side. That will be this, if edge_dist==i
      return ( edge_dist == n1->i-1 );
    }
    else return ( edge_dist < edge_dist2 );
  }
  else return ( (n1->j-n1->i) < (n2->j-n2->i) );
}


/**
 * Initialize the front
 */
void
Node::NewFront()
{
  for (int i=0; i<matrix_size+1; i++) Node::front.push_back(std::vector<bool>(matrix_size+1, false));
  
  Node::front_energy=EvalEnergy(Node::sequence.substr(0,Node::front_structure.size()),Node::front_structure);  

}

/**
 * Print matrix coordinates.
 */
inline std::string
Node::Print(Node* n)
{
  return ("(" + Str(n->i) + "," + Str(n->j) + ")");
}

/**
 * Print a vector of Nodes.
 */
std::string
Node::Print(std::vector<Node*> nodes)
{
  std::string s;
  for (unsigned int i=0; i<nodes.size(); i++) s += Print(nodes[i]);

  return (s);
}

/**
 * Print the front as a n+1xn matrix with 1's indicating squares that
 * are covered by the front.
 */
std::string
Node::PrintFront()
{
  std::string s;
  for (size_t i=0; i<front.size(); i++) {
    std::string nix(i, ' ');

    if ( i > 0 )s += nix;
    for (size_t j=i; j<front[i].size()-1; j++)
      s += Str(front[i][j]);
    s += "\n";
  }

  return (s);
}

/**
 * Removes all elements from Node::extrema that do not extend the
 * energy front.
*/

void
Node::PruneExtrema()
{
  std::vector<Node*> new_extrema=std::vector<Node*>();
  for(size_t i=0;i<extrema.size();i++){
    if(extrema[i]->FrontIncreases()){
         new_extrema.push_back(extrema[i]);
    }
  }
  extrema=new_extrema;
  
}

bool Node::Contains(Node* n) { return ((i <= n->i) && (j >= n->j));}


 void Node::ConstructFrontResolution(std::string & structure,std::vector<std::pair<int,int> > & basepairs, std::pair<int,int> bp_to_add){
  bool debug_pt=false;
  if(debug_pt) {
    Cout("ConstructFrontResolution"); //std::cout<<std::endl;
    std::cout<<"Basepairs "+PrintBasePairList(basepairs);
    std::string old_pair_table=PrintPairTable();
    Cout("old_pair_table:"+old_pair_table);
  }
  for(std::vector<std::pair<int,int> >::iterator it2=basepairs.begin();it2!=basepairs.end();it2++) {

    if(debug_pt) Cout("Conflict between "+PrintBasePair(*it2)+" and "+PrintBasePair(bp_to_add)+"?\n");
    if( BasePairConflict()(*it2, bp_to_add)){
      if(debug_pt) Cout("Remove "+PrintBasePair(*it2)+"\n");//std::cout<<std::endl;
      int first=it2->first;
      int second=it2->second;
      pair_table[first]=0;
      pair_table[second]=0;
      structure[first-1]='.';
      structure[second-1]='.';
      // bp_removed.push_back(*it2);
    }
  }


  std::vector<std::pair<int,int> >::iterator it=remove_if(basepairs.begin(),basepairs.end(),std::bind2nd(BasePairConflict(),bp_to_add));
  basepairs.erase(it,basepairs.end());
  if(debug_pt) Cout(PrintBasePairList(basepairs)+"\n");
  if(debug_pt) Cout(PrintPairTable());
 
  //add the basepair in stack, which now does not conflict with the front any more.
  basepairs.push_back(bp_to_add);
  if(debug_pt) Cout("Add "+PrintBasePair(bp_to_add)+"\n");
  int first=bp_to_add.first;
  int second=bp_to_add.second;
  pair_table[first]=second;
  pair_table[second]=first;
  structure[first-1]='(';
  structure[second-1]=')';
  if(debug_pt) Cout(PrintPairTable());
  if(debug_pt) Cout("new_pair_table:"+PrintPairTable());
}


void Node::BacktrackFront(std::pair<double,std::string> & ret,std::string node_substructure)
{
 
  bool debug_pt=false;//true;//false;
  if(debug_pt)  Cout("BacktrackFront with PairTable "+PrintPairTable());
  if(debug_pt)  Cout("node_substructure\n "+node_substructure+"\n");
  std::string front_structure=Node::front_structure;
  //Cout("have front_structure "+front_structure+" sequence "+sequence);
  //initialized later when front_structure with resolved conflicts is evaluated
  double front_energy;
  //create a substructure s
  //1.For a given Node n, add as much of the suboptimal structure s(n) of n 
  //2.Add more of it, if it improves E(s)
  //3.s=s+s(n)  
    std::vector<std::pair<int,int> >basepairs=MakeBasePairList1(front_structure);
    std::vector<std::pair<int,int> > bp_node=MakeBasePairList1(node_substructure);
    //take maximal possible substructure of node_substructure that does not conflict with front_structure
    std::vector<std::pair<int,int> >bp_node_conflict=std::vector<std::pair<int,int> >();
    for(size_t i=0;i<bp_node.size();i++){
      //do nothing if the base pair is already part of basepairs
       if(find(basepairs.begin(),basepairs.end(),bp_node[i])!=basepairs.end()) continue;
      if(!Conflict(basepairs,bp_node[i])) {
        //if (debug) {
        if(debug_pt) Cout("added bp "+PrintBasePair(bp_node[i])+"\n");//}
        basepairs.push_back(bp_node[i]);
        int first=bp_node[i].first;
        int second=bp_node[i].second;
        pair_table[first]=second;
        pair_table[second]=first;
        //manually change front_structure to circumvent a call to BasePairListToStructure1
        front_structure[first-1]='(';
        front_structure[second-1]=')';
      }
      else bp_node_conflict.push_back(bp_node[i]);
    }
    if(debug_pt) Cout("PairTable after adding nc bp \n"+PrintPairTable());
  

    //now basepairs is the maximal substructure without conflicts.
    //current substructure for the front
    //front_structure=BasePairListToStructure1(Node::sequence.size(),basepairs);
    front_energy=FastEvaluate();
    if(debug_pt) Cout("Front energy after adding non-conflicting base pairs is "+Str(front_energy)+" with PairTable \n"+PrintPairTable());
    if(verbose >= 3) Cout("front after adding non conflicting base pairs\n"+front_structure+":"+Str(front_energy)+"\n");

    //now add the baspairs in added_node that conflict with the front (and remove those basepairs in the front that conflict with them).
    //front after those basepairs in conflict with a basepair in added_note have been removed
    std::vector<std::pair<int,int> > new_basepairs= std::vector<std::pair<int,int> >();
    int min_stacksize=1;
    std::vector<std::vector<std::pair<int,int> > > bp_node_conflict_stacks=std::vector<std::vector<std::pair<int,int> > >();
    ConformationToStacks(bp_node_conflict_stacks,bp_node_conflict,min_stacksize);
    if(bp_node_conflict_stacks.size()>0){
      std::vector<std::pair<int,int> > bp_removed;
      for(size_t i=0;i<bp_node_conflict_stacks.size();i++){
        new_basepairs=basepairs;
        std::vector<std::pair<int,int> > stack=bp_node_conflict_stacks[i];
        //remove the basepairs in the stack that are already contained in the front
        for(int j=stack.size()-1;j>=0;j--) {
          if(ConformationHasPair(basepairs,stack[j])) stack.erase(stack.begin()+j);
        }
       
        std::vector<std::pair<double,std::string> > front_resolutions= std::vector<std::pair<double,std::string> >();
        //std::vector<std::pair<double,std::vector<std::pair<int,int> > > > front_resolutions= std::vector<std::pair<double,std::vector<std::pair<int,int> > > >();

        //front_resolutions.push_back(make_pair(front_energy,basepairs));
        front_resolutions.push_back(make_pair(front_energy,front_structure));//basepairs));     

        //ensure that the base pair in front at first_idx_stack does not come from the node that stack is constructued from.
        std::string resolved_structure=front_structure;
        for(size_t j=0;j<stack.size();j++){
          //last_insertion=j;
          //double resolution_energy=
          ConstructFrontResolution(resolved_structure,new_basepairs,stack[j]);
          front_resolutions.push_back(make_pair(FastEvaluate(),resolved_structure));
          //front_resolutions.push_back(make_pair(resolution_energy,BasePairListToStructure1(Node::transcribed,new_basepairs)));
        } 
        //revert changes in pair_table: remove stack and add those elts of basepairs back that were removed because they conflicted with elements in stack
        MakePairTable(const_cast<char*>(front_structure.substr(0,Node::transcribed).c_str()));        
        if(debug_pt) std::cout<<"Pairtable after reverting changes of outside in "+PrintPairTable();
        //add stack back to front as much as possible
        new_basepairs=basepairs;
        resolved_structure=front_structure;
        
        for(int j=stack.size()-1;j>=0;j--){
          ConstructFrontResolution(resolved_structure,new_basepairs,stack[j]);
          front_resolutions.push_back(make_pair(FastEvaluate(),resolved_structure));
        }

        //now accept the best solution in front_resolutions  and update structure, energy and bp
        front_energy=INF;
        front_structure=std::string(Node::transcribed,'.');
        for(std::vector<std::pair<double,std::string> >::iterator it=front_resolutions.begin();it!=front_resolutions.end();it++){
          if(it->first<front_energy){
            front_energy=it->first;
            front_structure=it->second;  
          }
        }
        basepairs=MakeBasePairList1(front_structure);         
        //sync pair_table with latest front unless this was the last iteration or the latest front is already in pair_table
        if(i<bp_node_conflict_stacks.size()-1){
          MakePairTable(const_cast<char*>(front_structure.substr(0,Node::transcribed).c_str()));
          if(debug_pt) std::cout<<"Pairtable after syncing with front resolution "+PrintPairTable();
          if(debug_pt) std::cout<<"Basepairs "+PrintBasePairList(basepairs);
      }
      //repeat for all conflict stacks
      }
    }
    if(debug_pt) Cout("BacktrackFront returns "+front_structure+" "+Str(front_energy)+"\n");
    ret= make_pair(front_energy,front_structure);
}



std::string Node::BacktrackNode(Node* n){
  if(n->IsMfE()) {
    Node::CalculateMfe();
    return Node::mfe_structure;
  }
  else {
    char *s = backtrack_fold_from_pair(const_cast<char*>(Node::sequence.c_str()),
                                       n->i, n->j);
    std::string ss(s);
    free(s);
    return (ss.substr(0,Node::transcribed));
  }
}

void Node::CalculateMfe(){
 
  char *sequence  = new char[Node::matrix_size+1];
  char *structure = new char[Node::matrix_size+1];
  strcpy(sequence, Node::sequence.c_str());
  //  if(fold_constrained)  strcpy(structure,Node::constraint_string.c_str()); 
  Node::mfe = fold(sequence, structure);
  Node::mfe_structure=std::string(structure);
  // clean up memory
  delete[] sequence;
  delete[] structure;
}

bool Node::IsMfE(){
  if(this->i==1 && this->j==(int)Node::sequence.size()) return true;
  else return false;

}



void Node::Transcribe(){ 
  Node::transcribed++;
  if(Node::transcribed==Node::matrix_size ) Node::SetAllEligible();
  if(verbose>=2)  Cout("Transcribed to "+Str(Node::transcribed)+"\n");
  S[0]=Node::transcribed;
  S1[0]=Node::transcribed;
  pair_table[0]=Node::transcribed;
  Node::front_structure+=".";
  MakePairTableFromFrontStructure();
  Node::front_energy=Evaluate(Node::front_structure);
}

bool Node::IsTranscribed(){
  return j<=Node::transcribed;

}

/**
TO DO:CAN IT BE THAT STRUCTURES SMALLR THAN MFE MAY STILL CONTRIBUTE TO FRONT WHILE THE MFE CANT
*/

bool Node::StopCriterionReached(){
  if(front_energy-MIN_ENERGY_DIFF<mfe){
    if(verbose>=1) Cout("StopCriterion reached: MfE is part of the front\n");
    return true;
  }
  return false;
}

double Node::MinimalExtensionCost(){
   if(extension_cost.size()==0) return 0.0; 
   else return *min_element(Node::extension_cost.begin(),Node::extension_cost.end());
}

void Node::RaiseEnergyBarrier(){
   Node::energy_barrier += 1.0;
   double minimial_extension_cost=MinimalExtensionCost();
   if(Node::energy_barrier < minimial_extension_cost && minimial_extension_cost>=0.0) Node::energy_barrier = ceil(minimial_extension_cost);
   Node * n;
   int count=0;
   for(size_t i=0;i<extrema.size();i++){
     n=extrema[i];
     if(!n->is_included && !n->is_eligible  &&  extension_cost[i]<=Node::energy_barrier){
         n->is_eligible=true;
         count++;
     }
   }
}


void Node::IncreaseTime(double inc){
  double oldtime=Node::time;
  Node::time+=inc;
  if(verbose>=1)  std::cout<<"Increase time from "<<oldtime<<" to "<<Node::time<<std::endl;
}

void Node::SetEnergyBarrier(){
  //time in seconds it takes to transcribe one base
 double t_transcribe=1.0/Node::transcription_rate;
 //number of bases transribed since start of algorithm
 double transcribed=Node::transcribed-Node::initially_transcribed;
 //time when next base is transcribed
 double new_time=(transcribed+1)*t_transcribe;
 double delta_time=new_time-Node::time;
 /*
 Cout("t_transcribe :"+Str(t_transcribe ));
 Cout("transcribed :"+Str(transcribed ));
 Cout("new_time :"+Str(new_time ));
 std::cout<<"SetEnergyBarrier for delta time "<<delta_time<<std::endl;
 */
 Node::SetEnergyBarrier(delta_time);
}

void Node::SetEnergyBarrier(double delta_time){
  Node::energy_barrier=Node::TimeToDeltaG(delta_time);
  //Cout("delta_time: "+Str(delta_time));
  //Cout("energy barrier: "+Str(Node::energy_barrier)+"\n");
}


void Node::AdvanceTimeToNextTranscription(){
  double t_transcribe=(double)1.0/(double)Node::transcription_rate;
    double transcribed=Node::transcribed-Node::initially_transcribed;
    double new_time=(transcribed+1)*t_transcribe;
    double delta_time=new_time-Node::time;
    /*
    Cout("t_transcribe: "+Str(t_transcribe));
    Cout("transcribed: "+Str(transcribed));
    Cout("new_time: "+Str(new_time));
    Cout("delta_time: "+Str(delta_time));
    */
    if(delta_time<EPSILON) delta_time+=t_transcribe;
    //    Cout("AdvanceTimeToNextTranscription:"+Str(delta_time));
    Node::IncreaseTime(delta_time);
    Node::Transcribe();
}


double Node::Evaluate(std::string struc){
  struc=struc.substr(0,Node::transcribed);
  std::string sequence=Node::sequence.substr(0,Node::transcribed);
  //Cout("Evaluate:\n"+sequence+"\n"+struc+"\n");
  return EvalEnergy(sequence,struc);
}


double Node::FastEvaluate(){
  std::string sequence=Node::sequence.substr(0,Node::transcribed);
  //  Cout("Evaluate:\n"+sequence+"\n"+struc+"\n");
  return FastEvalEnergy(sequence);

}
void Node::CleanUpNode(void)
{
  // free allocated memory of extrema array
  for (std::vector<Node*>::iterator it = Node::extrema.begin();
       it !=  Node::extrema.end();
       it++)
    delete *it;

  delete [] pair_table;
  delete [] S;
  delete [] S1;
}


bool Node::IsEligible(){
  return is_eligible;
}


bool Node::IsIncluded(){
  return is_included;
}

void Node::SetIncluded(bool val){
  is_included=val;
}


void Node::SetAllEligible(){
  Node::ineligible_count=0;
  for(size_t i=0;i<extrema.size();i++){
    Node * n=extrema[i];
    if(!n->is_included) n->is_eligible=true;
    else  Node::ineligible_count++;
  }
}

bool Node::HaveOverlap(const std::vector<std::pair<int,int> >& bp_removed,const std::pair<int,int> & extent){

 for(size_t i=0;i<bp_removed.size();i++){
   if(extent.second>=bp_removed[i].first && extent.first<=bp_removed[i].second ) return true;
 }
 return false;

}

void Node::Reeligify(std::string old_structure,std::string new_structure){
  int count=0;
  std::vector<std::pair<int,int> >bp_old=MakeBasePairList1(old_structure);
  std::vector<std::pair<int,int> >bp_new=MakeBasePairList1(new_structure);

  std::vector<std::pair<int,int> > bp_removed=std::vector<std::pair<int,int> > ();
  bp_removed.push_back(std::pair<int,int>(this->i,this->j));
  for(size_t i=0;i<bp_old.size();i++){
    if(find(bp_new.begin(),bp_new.end(),bp_old[i])==bp_new.end()) bp_removed.push_back(bp_old[i]);
  }
  Node * n;
  for(size_t i=0;i<extrema.size();i++){
    n=extrema[i];
    if(n->is_included || n->IsEligible()) continue;
    std::pair<int,int> extent(n->i,n->j);
    if(HaveOverlap(bp_removed,extent)) {
      n->is_eligible=true;
      count++;
    }
  }
}

void Node::SetIneligible(){
  is_eligible=false;
  Node::ineligible_count++;
}
 
void Node::MakePairTableFromFrontStructure(){
  //  std::string old_pair_table=PrintPairTable();
  MakePairTable(const_cast<char*>(Node::front_structure.substr(0,Node::transcribed).c_str()));
  //   pair_table[0]=Node::transcribed;
  //  std::cout<<"MakePairTableFromFrontStructure "+Node::front_structure<<std::endl;
  //std::cout<<"old:"+old_pair_table<<std::endl;
  //std::cout<<"new:"+PrintPairTable()<<std::endl;
}

bool Node::AllTranscribed(){
  return (Node::transcribed==Node::matrix_size);
}
