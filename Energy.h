/*
  Last changed Time-stamp: <2006-02-16 10:49:19 xtof>
  $Id: Energy.h,v 1.4 2007/04/16 11:20:10 Kinwalker Exp $
*/
#ifndef _ENERGY_H_
#define _ENERGY_H_

#include <string>
#include <vector>
#include "options.h"
#include "Util.h"


extern "C"{
  //#include "read_epars.h"
 #include "utils.h"
#include "energy_const.h"
#include "energy_par.h"
  extern void  read_parameter_file(const char *);


  //extern 
  //paramT * P;
}


double/*float*/
EvalEnergy(std::string sequence, std::string structure);
double/*float*/
FastEvalEnergy(std::string sequence);
void
InitializeEnergyModel(OptionS* OptS, std::string sequence);
bool
CanPair(char i, char j);
std::vector<int>
EncodeSequence(std::string sequence);
void
InitializeViennaRNA(std::string sequence,int dangle,int transcribed);
float
FullEnergyModel(std::string sequence, std::string structure);
float
NJEnergyModel(std::string sequence, std::string structure);
int
NJ_energy(int i, int j);
std::vector<std::pair<float,std::string> >
EnergyEvalStructureList(std::string sequence,
			std::vector<std::string> strList);
std::vector<int>
MakeLoopIndex(std::string structure);
std::vector<int>
MakePairTable(std::string structure);
std::vector<std::pair<int,int> >
MakeBasePairList (std::string structure);
std::vector<std::pair<int,int> >
MakeBasePairList1 (std::string structure);

//extern float
//(*EvalEnergy)(std::string sequence, std::string structure);

#endif

/* End of file */
