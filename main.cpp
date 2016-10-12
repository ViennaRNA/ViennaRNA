/*
  Last changed Time-stamp: <2007-04-02 15:19:44 caro>
  $Id: main.cpp,v 1.21 2007/11/04 14:59:54 Kinwalker Exp $
*/
#include "Node.h"
#include <ctime>

extern "C" {
#include "fold.h"

}
//./kinwalker --barrier_heuristic=S --transcribed=1 --transcription_rate=200 --dangle=0 -vvv < test/sv11.seq
// ./kinwalker --transcribed=1 --transcription_rate=20 (--dangle=0/2) < test/kinfold_test.seq
//   RNAeval -d2 -P /home/mescalin/mgeis/kinwalker/noTermAU.par
// valgrind --tool=memcheck --leak-check=full ./kinwalker --transcribed=1 --transcription_rate=200 < test/kinfold_test.seq
  char *structure=NULL, *cstruc=NULL;

int main (int argc, char *argv[]) {



  //Transcription Rate in base per second
  //Eukaryonten 10-20
  //Bacteriophagen 200 
  //Bacterien 20-80

  // Process commandline options  
  OptionS* OptS;
  OptS = decodeCML(argc, argv);
 
  clock_t start,finish;
  double time;
  start = clock();
 

  Node::ProcessOptions(OptS);
  InitializeEnergyModel(OptS,Node::sequence);
  if(OptS->init_structure) Node:: trajectory.push_back(Node::front_structure+" "+Str(Node::Evaluate(Node::front_structure))+" "+Str(Node::time)+" "+Str(0)+" "+Str(0)+" "+Str(Node::transcribed));
 
  Node::CalculateMfe();
  Cout(Node::sequence+"\n");
  Cout(Node::mfe_structure+" "+Str(Node::mfe)+"\n");
 

  // Find Extrema and create successively greater fronts with them.
  Node::CalculateMaxBarrier();

  Node::FindLocalExtrema();

  Node::NewFront();  


  //////////////////////////////////
  //ALGORITHM:
  //Repeat the following until the final structure is reached
  //1.Try to Extend the Front before the next base is transcribed
  //2.If front could not be extended:
  //a. Transcribe one base unless sequence already is entirely transcribed
  //b. If entirely transcirbed, increase the energy barrier
  //////////////////////////////////


  //while kinwalker has not finished yet
  while (!Node::StopCriterionReached() ) {
    if(Node::transcribed<Node::matrix_size) Node::SetEnergyBarrier();
    //try to extend the front given the transcribed base pairs and the energy barrier
    bool extended=Node::ExtendFront();
    if(!extended) {
      //Cout("Not extended with transcribed "+Str(Node::transcribed)+" and matrix size "+Str(Node::matrix_size)+"\n");
      if(Node::transcribed<Node::matrix_size) {
        Node::AdvanceTimeToNextTranscription();       
        if(Node::verbose>=1) Cout(Str(Node::transcribed)+" bases transcribed\n");
      }
       //if everything is already transcribed, increase the energy barrier.
       else Node::RaiseEnergyBarrier(); 
    }
  }
  Cout("TRAJECTORY\n");
  for(size_t i=0;i<Node::trajectory.size();i++){
    Cout(Node::trajectory[i]+"\n");
  }
 
  //if(Node::verbose>=1) 
  finish = clock();
  time = (double(finish)-double(start))/CLOCKS_PER_SEC;
  Cout("Kinwalker run time:"+Str(time)+" seconds \n");

  // clean up memory!!!
  Node::CleanUpNode();
  if(Node::print_front_trajectory) {
    Node::front_trajectory_ps+="end \n";
   std::ofstream outfile("front_trajectory.ps");
   outfile << Node::front_trajectory_ps;
   outfile.close();
  }
  return (EXIT_SUCCESS);
}

// End of file
  
