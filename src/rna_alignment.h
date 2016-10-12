#ifndef _RNA_ALIGNMENT_H
#define _RNA_ALIGNMENT_H

#include <iostream>
#include <map>
#include <string>
#include "alignment.h"
#include "forestali.h"
#include "rna_algebra.h"
#include "rnaforest.h"
#include "rnafuncs.h"
#include "alignment.t.cpp"



// TODO name ist irrefuehrend..
// RNA alignment forest, inherits alignment algorithms, 
// and adds how to work with the labels, and stucture names, pair table and structure info
class RNA_Alignment : public ForestAli<RNA_Alphabet,RNA_AlphaPair> {
private:
    std::string strname1_;
    std::string strname2_;

    // Implementations of virtual functions

    void makeRepLabel(size_type node, RNA_Alphabet a,  RNA_Alphabet b, int anchor) {
        lb_[node].a=a;
        lb_[node].b=b;
				// anchor not needed in resulting alignment structure for pairw. ali
    };

    void makeDelLabel(size_type node, RNA_Alphabet a) {
        lb_[node].a=a;
        lb_[node].b=ALPHA_GAP;
    };

    void makeInsLabel(size_type node, RNA_Alphabet b) {
        lb_[node].a=ALPHA_GAP;
        lb_[node].b=b;
    };

    void showLabel(std::ostream &s,RNA_AlphaPair p) const {
        s << "[" << p.a << p.b << "]";
    };


    void makePairTable(std::map<unsigned int,unsigned int> &pairs, bool first) const;

public:
    void setStructureNames(const std::string &s1,const std::string &s2);
    const std::string& getStructureNameX() const {
        return strname1_;
    };
    const std::string& getStructureNameY() const {
        return strname2_;
    };
    void getSequenceAlignments(std::string &s1, std::string &s2) const;
    void getStructureAlignment(std::string &s, bool first) const;
    void generateXML(std::ostream &s) const;
    void squigglePlot(const std::string &filename_suffix, const RNAFuncs::SquigglePlotOptions &options) const;
};

#endif
