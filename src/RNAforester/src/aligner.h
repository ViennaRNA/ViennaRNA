#ifndef ALIGNER_H
#define ALIGNER_H

#include "alignment.h"
#include "debug.h"
#include "treeedit.h"

#include "misc.h"
#include "progressive_align.h"
#include "rna_alignment.h"
#include "rnaforest.h"
#include "rnaforestsz.h"
#include "rnafuncs.h"
#include "rna_profile_alignment.h"
#include "rna_algebra.h" // TODO score rausholen in eigenes file

#include "alignment.t.cpp"
#include "treeedit.t.cpp"

#include <vector>
#include <sys/times.h>


// forward decls
class RNAForest;

// TODO diese freunde gehoeren eigentl. in printer
template <class L>
void makeDotFileAli(const Forest<L> &f, const Options &options) {
    if (options.has(Options::OutputAlignmentDotFormat)) {
        std::string filename;
        options.get(Options::OutputAlignmentDotFormat,filename,std::string("ali.dot"));
        std::ofstream s(filename.c_str());
        f.printDot(s);
    }
}

template <class L>
void makeDotFileInp(const Forest<L> &f, const Options &options, unsigned int count) {
    if (options.has(Options::MakeDotForInputTrees)) {
        std::ostringstream ss;
        std::ofstream os;
        ss << "input" << count << ".dot";
        os.open(ss.str().c_str());
        f.printDot(os);
        os.close();
    }
}

void editPairwise(std::vector<RNAForestSZ*> &inputListSZ, Score &score, Options &options, bool anchored);
void alignPairwiseSimple(std::vector<RNAForest*> &inputListPW, Score &score, Options &options, bool anchored);
void alignPairwise(std::vector<RNAForest*> &inputListPW, Score &score, const Options &options, bool anchored, RNAFuncs::AddXmlInfos &xmlInfos);
void computeSpaceTimeInfo(const RNAForest &f1, const RNAForest &f2, const Score &score, const Options &options);

void alignMultiple(std::vector<RNAProfileAlignment*> &inputListMult, Score &score, const Options &options, bool anchored, RNAFuncs::AddXmlInfos &xmlInfos);
void printCluster(std::pair<double, RNAProfileAlignment*> cluster,
                  unsigned int clusterNr, double minPairFreq, const Options & options);

inline double getAndPrintOptimum(Alignment<double,RNA_Alphabet,RNA_AlphaPair> * ali, const Options & options); 

#endif  /* _ALIGNER_H */
