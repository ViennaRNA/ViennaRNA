/*
  Last changed Time-stamp: <2006-10-13 22:34:37 xtof>
  $Id: Energy.cpp,v 1.12 2007/11/03 16:45:58 Kinwalker Exp $
*/
#include <cstdlib>
#include <iostream>
#include <string>
#include <cstring>
#include "Energy.h"


extern "C" {
#include "fold.h"
#include "fold_vars.h"
#include "utils.h"
#include "pair_mat.h"
}
extern short * S;
extern short * S1;
extern short * pair_table;

static float
(*EnergyModel)(std::string sequence, std::string structure) = NULL;




double
FastEvalEnergy(std::string sequence){
  return energy_of_struct_pt(const_cast<char*>(sequence.c_str()), pair_table, S, S1)/100.;
  //= energy_of_struct(sequence.c_str(), structure.c_str());
}
//operation: einfuege, wegnehmen
//jetzige gesamtenergue
//welche loops betroffen
//inline int  LoopEnergy(int n1, int n2, int type, int type_2,int si1, int sj1, int sp1, int sq1);


//   void  PartialPath::RemoveFromPairTableWithLog(int first,int second,int removeCatalogIdx){
//     //    cout<<"RemoveFromPairTableWithLog idx "+Str(removeCatalogIdx)<<endl;
//     currentPartialPathTrail.push_back(removeCatalogIdx);
//     RemoveFromPairTable(first,second);
//   }

//  void  PartialPath::AddToPairTableWithLog(int first,int second,int addCatalogIdx){
//    //   cout<<"AddToPairTableWithLog idx "+Str(addCatalogIdx)<<endl;
//     currentPartialPathTrail.push_back(addCatalogIdx+BP_ADD_CONST);
//     AddToPairTable(first,second);
//   }



//   void PartialPath::AddToPairTable(std::pair<int,int> bp){
//     // cout<<"AddToPairTable"+PrintBasePair(bp) <<std::endl;
//     pair_table[bp.first]=bp.second;
//     pair_table[bp.second]=bp.first;
//   }


//   void PartialPath::RemoveFromPairTable(std::pair<int,int> bp){
//     // cout<<"RemoveFromPairTable" +PrintBasePair(bp)<<std::endl;
//     pair_table[bp.first]=0;
//     pair_table[bp.second]=0;
//   }




double
EvalEnergy(std::string sequence, std::string structure)
{
  bool debug=false;
  if(debug){
  std::cout<<"EvalEnergy"<<std::endl;
  std::cout<<sequence<<std::endl;
  std::cout<<structure<<std::endl;
  }
  return (EnergyModel(sequence, structure));
}


void
InitializeEnergyModel(OptionS* OptS, std::string sequence)
{
  //hard code to ViennaRNA as engies are being minimized in lots of places
    InitializeViennaRNA(sequence,OptS->dangle,OptS->transcribed);
    EnergyModel = FullEnergyModel;

 //  if (OptS->Emodel == 'M') {
//     InitializeViennaRNA(sequence,OptS->dangle,OptS->transcribed);
//     EnergyModel = FullEnergyModel;
//   }
//   else
//     EnergyModel = NJEnergyModel;
}


bool
CanPair(char i, char j)
{
  return (pair[encode_char(i)][encode_char(j)] > 0);
}

/**
 * Numerically encodes the sequence
 */
std::vector<int>
EncodeSequence(std::string sequence)
{
  std::vector<int> encoded;

  for (std::string::const_iterator base = sequence.begin();
       base != sequence.end();
       base++)
    encoded.push_back(static_cast<int>(encode_char(*base)));

  return (encoded);
}

/**
 *
 */



void
InitializeViennaRNA(std::string sequence,int dangle,int transcribed)
{
  int length=sequence.size();
  S  = new short [length+2];
  S1  = new short [length+2];//NULL;//new short [length+2];
  pair_table =  new short [length+2];
  // pair_table = make_pair_table(const_cast<char*>(std::string(length,'.').c_str()));
  //(short *) MG_space(sizeof(short)*(length+1));
  MakePairTable(const_cast<char*>(std::string(length,'.').c_str()));
  make_pair_matrix();
  S[0] = S1[0] = pair_table[0]=transcribed;
  for (int i=0; i< length; i++) {
    S[i+1] = encode_char(sequence[i]);
    S1[i+1] = alias[S[i+1]];
  } 

  //  std::cout<<"Dangle "<<dangle<<std::endl;
  dangles = dangle;//2;
  //  no_closingGU=1; no effect 
  
  initialize_fold(length);
  update_fold_params();
  //  read_parameter_file("/home/mescalin/mgeis/kinwalker/noTermAU.par"); no effect 
}

/**
 * Calculates the energy of a nucleic acid using the energy model as
 * implemented in the Vienna RNA package
 */
float
FullEnergyModel(std::string sequence, std::string structure)
{
  // initialize_fold(sequence.length());
 
  float energy = energy_of_struct(sequence.c_str(), structure.c_str());
  // free_arrays();

  return (energy);
}

/**
 * Calculates the energy of a nucleic acid using the Nussinov-Jacobson Model
 */
float
NJEnergyModel(std::string sequence, std::string structure)
{
  int energy = 0;
  int length = structure.length();
  std::vector<int> sq   = EncodeSequence(sequence);
  std::vector<int> ptbl = MakePairTable(structure);

  for (int i=0; i<length; i++) {

    if ( ptbl[i] == -1 ) continue; // unpaird positions
    if ( ptbl[i] < i ) continue;   // positions with j<i
    energy += NJ_energy(sq[i], sq[ptbl[i]]);
  }

  return (static_cast<float>(energy));
}

/**
 * Base pair energies rules of the Nussinov-Jacobson Model
 */
int
NJ_energy(int i, int j)
{
  int A, U, G, C;
  A = encode_char('A');
  U = encode_char('U');
  G = encode_char('G');
  C = encode_char('C');

  if ( (i == G && j == C) || (i == C && j == G) ) return (3); /* GC CG */
  if ( (i == A && j == U) || (i == U && j == A) ) return (2); /* AU UA */
  if ( (i == G && j == U) || (i == U && j == G) ) return (1); /* GU UG */
  return (0); /* unknown base pair */
}

/**
 * Calculates the energies for a list of structures
 * [Structure] => [(Enegy, Structure)]
 */
std::vector<std::pair<float,std::string> >
EnergyEvalStructureList(std::string sequence, std::vector<std::string> strList)
{
  std::vector<std::pair<float,std::string> > tupelList;

  for (std::vector<std::string>::const_iterator structure = strList.begin();
       structure != strList.end();
       structure++) {

    float energy = FullEnergyModel(sequence, *structure);
    tupelList.push_back(make_pair(energy, *structure));
  }

  return (tupelList);
}

/**
 * Number loops and assign each sequence position its loop-number
 */
std::vector<int>
MakeLoopIndex(std::string structure)
{
  int hx, l, nl;
  int length;
  std::vector<int> stack;
  std::vector<int> loop;

  length = structure.size();
  stack  = std::vector<int>(length+1, -1);
  loop   = std::vector<int>(length+2, -1);

  hx = l = nl = 0;
  for (int i=0; i<length; i++) {
    if ( structure[i] == '(' ) {
      nl++; l = nl;
      stack[hx++] = i;
    }
    loop[i] = l;
    if ( structure[i] ==')' ) {
      --hx;
      if ( hx > 0 )
	l = loop[stack[hx-1]];  // index of enclosing loop
      else l = 0;               // external loop has index 0
      if ( hx < 0 ) {
     Cout("MakeLoopIndex MakePairTable(): "+structure+" unbalanced brackets");
     exit(EXIT_FAILURE);
       }
      // Fatal("MakePairTable(): %s\n%s\n", "unbalanced brackets", structure.c_str());
    }
  }

  return (loop);
}

#if 0

/**
 * let table[i]=j if (i,j) is pair, -1 if i unpaired indices start at
 * 0 in this version!
 */
std::vector<int>
MakePairTable(std::string structure)
{
  short* pt = NULL;
  std::vector<int> tbl;

  pt = make_pair_table(structure.c_str());
  for(int i = 1; i<=pt[0]; i++) {
    int value = (pt[i]==0) ? -1 : static_cast<int>(pt[i]-1);
    tbl.push_back(value);
  }
  free(pt);
  return (tbl);
}

#else

/**
 * let table[i]=j if (i,j) is pair, -1 if i unpaired indices start at
 * 0 in this version!
 */
std::vector<int>
MakePairTable(std::string structure)
{
  //Cout("MakePairTable(): "+structure+"\n");
  int length = structure.size();
  std::vector<int> stk = std::vector<int>();
  std::vector<int> tbl(length, -1);

  for (int i=0,j=-1; i<length; i++) {
    if ( structure[i] == '(' )
      stk.push_back(i);
    else if ( structure[i] == ')' ) {
	j = stk[stk.size()-1];
	stk.pop_back();
	if ( j < 0 ) {
           Cout("MakePairTable(): "+structure+" unbalanced brackets");
          std::string s=std::string();
          for(int r=0;r<=tbl[0];r++){
            s+=Str(tbl[r])+" ";
          }
           Cout(s);
           exit(EXIT_FAILURE);
    }
	//	if ( j < 0 ) Fatal("MakePairTable(): %s\n%s\n",
	//	   "unbalanced brackets",
	//	   structure.c_str());
      
	tbl[i] = j;
	tbl[j] = i;
    }
  }

  if (! stk.empty() ) {
     Cout("MakePairTable(): "+structure+" unbalanced brackets stk not empty");
     exit(EXIT_FAILURE);
  }
    //Fatal("MakePairTable(): %s\n%s\n", "too few closing brackets", structure.c_str());

  return (tbl);
}

#endif

std::vector<std::pair<int,int> >
MakeBasePairList (std::string structure)
{
  std::vector<int> pTbl = MakePairTable(structure); 
  std::vector<std::pair<int,int> > pLst = std::vector<std::pair<int,int> > ();

  int i = 0;
  for (std::vector<int>::iterator j = pTbl.begin();
       j != pTbl.end();
       j++, i++) {

    // skip unpaired positions and entries with pTbl[i] < i
    if (*j == -1 || *j < i) continue;

    pLst.push_back(std::make_pair(i, *j));
  }

  return (pLst);
}

std::vector<std::pair<int,int> >
MakeBasePairList1 (std::string structure)
{
  std::vector<int> pTbl = MakePairTable(structure); 
  std::vector<std::pair<int,int> > pLst = std::vector<std::pair<int,int> > ();

  int i = 0;
  for (std::vector<int>::iterator j = pTbl.begin();
       j != pTbl.end();
       j++, i++) {

    // skip unpaired positions and entries with pTbl[i] < i
    if (*j == -1 || *j < i) continue;

    pLst.push_back(std::make_pair(i+1, *j+1));
  }

  return (pLst);
}
  


// End of file
