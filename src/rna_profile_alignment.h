#ifndef _RNA_PROFILE_ALIGNMENT_H
#define _RNA_PROFILE_ALIGNMENT_H

#include <vector>
#include <iostream>
#include <map>
#include "alignment.h"
#include "matrix.h"
#include "forestali.h"
#include "rna_algebra.h"


// TODO name ist irrefuehrend
// TODO squiggle plot und print etc. auslagern in print klasse!
// profile alignment forest 

class RNAProfileAlignment
            : public ForestAli<RNA_Alphabet_Profile,RNA_Alphabet_Profile> {
public:
    struct SquigglePlotOptions {
        bool hideBaseNumbers;
        unsigned int baseNumInterval;
        bool greyColors;
        bool mostLikelySequence;
        double minPairProb;
    };

    struct BaseProbs {
        double a;
        double c;
        double g;
        double u;
        double gap;
        double base;
    };

private:
    std::string name_;
    unsigned int numStructures_;
    unsigned int numStructuresX_;
    unsigned int numStructuresY_;
    std::vector<std::string> strNames_;
    bool hasSequence_;

    void makeRepLabel(size_type node, RNA_Alphabet_Profile a,  RNA_Alphabet_Profile b, int anchor);
    void makeDelLabel(size_type node, RNA_Alphabet_Profile a);
    void makeInsLabel(size_type node, RNA_Alphabet_Profile b);

    void showLabel(std::ostream &s,RNA_Alphabet_Profile p) const;
    void makeLabel(RNA_Alphabet_Profile &a,char c);

    /** Build forest from structure, sequence pair. */
    void buildForest(const std::string &baseStr, const std::string &viennaStr, bool use_bp_prob=false);
    void buildForestAnchored(const std::string &baseStr, const std::string &viennaStr, bool use_bp_prob=false);

    void getStructureAlignmentFromCSF(std::string &s, std::vector<double> &pairprop, double t,size_type i, size_type j) const;
    double bestPairs(size_type node) const;
    void drawBaseCircles(int device_id,const BaseProbs &bp,double center_x,double center_y) const;
    double getMlBaseFreq(const BaseProbs &baseprobs) const;
    char getMlBase(const BaseProbs &bp) const;
    void getSeqAli(std::string &seq,unsigned int row,unsigned int i,unsigned int j) const;
    void getStructAli(std::string &str,unsigned int row) const;
    void makePairTable(std::map<unsigned int,unsigned int> &pairs, unsigned int row) const;
    void filterConsensus(std::string &structure, std::vector<double> &pairprob, std::vector<BaseProbs> &baseprobs, double minFreq) const;
    void addStrName(const std::string &strName) {
        strNames_.push_back(strName);
    };

    inline bool isBase(size_type node) const {
        if (label(node).p[ALPHA_PRO_BASE] > 0)
            return true;
        else
            return false;
    };

    inline bool isPair(size_type node) const {
        return !isBase(node);
    };

public:
		// normal prof. ali
    RNAProfileAlignment(const std::string &baseStr, const std::string &viennaStr, const std::string &name, bool anchored);
		// for predict profile 
    RNAProfileAlignment(const std::string &baseStr, const std::string &constraint, const std::string &name, bool anchored, double t);
    RNAProfileAlignment(const std::string &filename, const bool anchored);
#ifdef HAVE_LIBRNA
		RNAProfileAlignment(const std::string &baseStr, const std::string &name, bool anchored, const std::string &constraint, double t);
#endif
		virtual ~RNAProfileAlignment(){ delete[] anchors_; };

    void printSeqAli() const;
    std::vector<std::pair<std::string,std::string> > getSeqAli() const;
    std::vector<std::string> getStrAli() const;
    std::string getConsSeq() const;
    std::string getConsStr(double minPairProb) const;
    std::vector<double> getBaseProb() const;
    void getPairProb(double &minPairProb, std::vector<std::pair<std::pair<int,int>,double> > &pairprobs);
    void printStrAli() const;
    void printConsensus(double minPairProb) const;
    void printFastaAli(bool noStructure=false) const;

    void squigglePlot(const std::string &filename, SquigglePlotOptions &options) const;
    void getSequenceAlignment(std::vector<BaseProbs> &baseprob) const;
    void getStructureAlignment(double t, std::string &s, std::vector<double> &pairprob) const;
    std::string getName() const {
        return name_;
    };
    void setName(const std::string &name) {
        name_ = name;
    };
    unsigned int getNumStructures() const {
        return numStructures_;
    };
    const std::vector<std::string>& getStrNames() const {
        return strNames_;
    };
    void addStrNames(const std::vector<std::string>& strNames);
    void save(const std::string &filename);

    RNAProfileAlignment(unsigned int numStructuresX,unsigned int numStructuresY) :
            ForestAli<RNA_Alphabet_Profile,RNA_Alphabet_Profile>(),
            numStructures_(numStructuresX + numStructuresY),
            numStructuresX_(numStructuresX),
            numStructuresY_(numStructuresY) {
						};
};

#endif
