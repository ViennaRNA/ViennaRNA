/*
  Last changed Time-stamp: <2007-07-10 14:09:36 xtof>
  $Id: Node.h,v 1.20 2007/11/03 16:45:58 Kinwalker Exp $
*/
#ifndef _NODE_H_
#define _NODE_H_

/*
#if __cplusplus
#  define BEGIN_C_DECLS extern "C" {
#  define END_C_DECLS   }
#else
#  define BEGIN_C_DECLS
#  define END_C_DECLS
#endif
*/


#ifdef __cplusplus


#include "options.h"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <cstring>
#include <vector>


extern "C" {
  #include "findpath.h"  
  #include "fold.h"
  #include "energy_const.h"
  #include "utils.h"
  #include "fold_vars.h"
}
#include "Energy.h"
#include "MorganHiggs.h"
#include "Util.h"





/**
 * Each Node is unique as matrix positions are unique. All
 * calculations of front etc. on the max matching matrix are included
 * in the Node class. Usually as static members.
*/
class Node {
 public:
static  std::vector<std::pair<double,std::string> > path;


  /*
set of time and structure pairs
  */

  static std::vector<std::string> trajectory;
  static int ineligible_count;
  
  //  static std::vector<std::pair<double,std::string> > trajectory;
  int i;
  int j;
  double free_energy;
  double saddle_barrier;
  std::string structure;
  bool is_eligible;
  bool is_included;

  //double rel_free_energy;
  //int pairs;
  //static bool done_once;
  static OptionS* OptS;
  static int verbose;
  static std::string grouping;
  static std::string sequence; 
  static double front_energy;
  static std::string front_structure;
  static std::string mfe_structure;
  static double mfe;
  static int lookahead;
  static int matrix_size;
  static double energy_barrier;
  static double max_barrier;
  static std::vector<std::vector<bool> > front;
  static std::vector<Node*> front_extrema;
  static std::vector<Node*> extrema;
  static std::vector<double> extension_cost;
  static double time;
  static double transcription_rate;
  static int transcribed;
  static int initially_transcribed;
  static bool transcribed_front_maximal;
  static int interrupt_trajectory;
  static std::string front_trajectory_ps;
  static bool print_front_trajectory;
  static std::string constraint_string;
  static int windowsize;
  /**
   * Constructor
   */
  Node(int i,int j, double energy);

  bool IsWithinWindow();

  static bool AllTranscribed();
  void AddToFront(std::string new_front_structure,double barrier);
  void AddToFront(std::string new_front_structure,double new_front_energy,double barrier);
  static bool FindExtremum();
  static bool AddsMoreToFront(Node* n1, Node* n2);
  void AddToFront();
  static void BacktrackFront(std::pair<double,std::string>& ret,std::string node_substructure);
  static std::string BacktrackNode(Node* n);
  static void CalculateMfe();
  bool Contains(Node* n);
  static void CalculateMaxBarrier();
  static bool ExtendFront();
  static void FindLocalExtrema();
  static std::string FoldingPathToString(std::vector<std::pair<double,std::string> > path);

  static bool StructureIsFrontExtension( const std::pair<double,std::string> & struc);
  std::string AddTerminalBasePairsToFrontStructure();
   
  static std::string PrintLocalExtrema();
  bool IsMfE();
  int FrontIncrease();
  bool FrontIncreases();
  bool IsTranscribed();
  static void ProcessOptions(OptionS* OptS);
  static void RaiseEnergyBarrier();
  static bool StopCriterionReached();
  static void Transcribe();
  bool TrulyContainsFront();
  static inline bool LessThan(Node* n1, Node* n2);
  static void NewFront();
  static inline std::string Print(Node* node);
  static std::string Print(std::vector<Node*> nodes);
  static std::string PrintFront();
  static void PruneExtrema();
  double EnergyBarrier(std::string structure);
  void CalculateSaddleBarrier(std::string structure);
  double EnergyBarrierFromStructure(std::string structure);
  static double MinimalExtensionCost();
  static double TimePassedFolding(double delta_G);
  static double DeltaGToTimeUnits(double delta_G);
  static double TimeUnitsToDeltaG(double units);
  static void IncreaseTime(double inc);
  static void SetEnergyBarrier();
  static void SetEnergyBarrier(double delta_time);
  static void AdvanceTimeToNextTranscription();
  static double TimeToDeltaG(double time);
  static double Evaluate(std::string structure);
  static double FastEvaluate();
 
  static void ConstructFrontResolution(std::string & structure,std::vector<std::pair<int,int> > & basepairs, std::pair<int,int> bp_to_add);

  static void CleanUpNode(void);

  static void MakePairTableFromFrontStructure();
  bool IsEligible();
  void SetIneligible();
  bool IsTemporarilyIneligible();
  void SetTemporarilyIneligible();
  bool IsIncluded();
  void SetIncluded(bool val);

  static void SetAllEligible();
  void Reeligify(std::string old_structure,std::string new_structure);
  static  bool HaveOverlap(const std::vector<std::pair<int,int> >& bp_removed,const std::pair<int,int> & extent);
  static void  UpdateEligibility();

    static void CalculateFoldingPath(Node * extremum,std::string node_substructure);

    static void GetSaddleFromPath(std::pair<double,std::string> & saddle,std::pair<double,std::string> & final_structure);
    static void GetMorganHiggsPath(std::string target );
    static int FindPathSaddle(int min_idx);
    static int FindPathMinimum(int last_idx_within_barrier);
    static std::string CombineFrontAndNode(std::string node_substructure);


    static void GetMorganHiggsStudlaPath(std::string target);

    static void CalculatePathAndSaddle(Node* extremum,std::string integrated_structure,std::pair<double,std::string> & saddle,std::pair<double,std::string> & final_structure);

};


#endif /* __cplusplus */
#endif  /* _NODE_H_ */
