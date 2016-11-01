#include "graphtypes.h"
#include "matrix.h"
#include "rnaforester_options.h"
#include "rna_profile_alignment.h"
#include "wmatch.h"

/* ****************************************** */
/*          Definitions and typedefs          */
/* ****************************************** */

typedef std::map<long,RNAProfileAlignment*> RNAProfileAliMapType;
typedef std::pair<long,RNAProfileAlignment*> RNAProfileAliKeyPairType;

/* ****************************************** */
/*            Function prototypes             */
/* ****************************************** */

void progressiveAlign(std::vector<RNAProfileAlignment*> &inputList, 
											std::vector<std::pair<double,RNAProfileAlignment*> > &resultList, const Score &score, const Options &options, bool anchoring);
Graph makePairsGraph(const RNAProfileAliMapType &inputListProfile, const Algebra<double,RNA_Alphabet_Profile> *alg, const Matrix<double> *score_mtrx, double threshold);
