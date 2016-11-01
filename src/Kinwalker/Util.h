#include <string>
#include <cstring>
#include <vector>
#include <sstream>
#include <cstdarg>
#include <iostream>
#include <cstdlib>
#include <algorithm>

void Cout(std::string s);

//void Fatal(char *fmt, ...);

std::string Str(double x);
std::string Str(float x);
std::string Str(int x);
std::string Str(unsigned int x);

 std::vector<std::pair<int,int> >
PairTableToBasePairList(short int * pair_table);

//void
//PrintBasePairList(std::vector<std::pair<int,int> > bpl);
std::string
BasePairListToStructure(int lenght, std::vector<std::pair<int,int> > pl);
std::string
BasePairListToStructure1(int length, const std::vector<std::pair<int,int> > & pl);
//std::string
//BasePairListToStructure1(int length, const std::vector<std::pair<int,int> > & pl);


std::string PrintBasePair(std::pair<int,int> bp);
std::string PrintBasePairList(std::vector<std::pair<int,int> > list);


bool IntroducesPseudoKnot(const std::vector<std::pair<int,int> >& node, const std::pair<int,int>& p1);

bool ConformationHasPair(const std::vector<std::pair<int,int> >& node,const std::pair<int,int> & p);

bool Conflict(const std::vector<std::pair<int,int> >& node, const std::pair<int,int>& p1);

bool Conflict(const std::vector<std::pair<int,int> >& node1, const std::vector<std::pair<int,int> >& node2);

//struct BasePairConflict :public std::binary_function<std::pair<int,int>,std::pair<int,int>,bool>;

//std::vector<std::vector<std::pair<int,int> > > ConformationToStacks(std::vector<std::pair<int,int> > node,int stacksize);

void
ConformationToStacks(std::vector<std::vector<std::pair<int,int> > > & stacks, std::vector<std::pair<int,int> > node,int stacksize);

std::string PrintPairTable();

void MakePairTable(const char *structure);

std::string PSFrontPlot(std::string sequence,std::vector<std::pair<int,int> > extrema);
