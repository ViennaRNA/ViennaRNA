#include <algorithm>
#include <string>

#ifndef WIN32
#include "config.h"
#endif

#include "misc.h"
#include "rna_alignment.h"
#include "rnafuncs.h"
#include "utils.h"

/* ****************************************** */
/*            Private functions               */
/* ****************************************** */

void RNA_Alignment::makePairTable(std::map<unsigned int,unsigned int> &pairs, bool first) const {
    std::pair<int,int> *baseIndex;   // stores for each node in the alignment the index position of the left and right base (or -1)
    RNA_Alphabet c;

    baseIndex = new std::pair<int,int>[size_];

    // initialize pairs, all gaps
    for (int i=size()-1; i>=0; i--) {
        baseIndex[i].first=-1;
        baseIndex[i].second=-1;
    }

    for (int i=size()-1; i>=0; i--) {
        // first or second component
        if (first)
            c=lb_[i].a;
        else
            c=lb_[i].b;

        if (isLeaf(i)) {
            if (c!=ALPHA_GAP) {
                baseIndex[i].first=i;
                baseIndex[i].second=i;
            }
        } else {
            // internal node
            // leftmost and rightmost base
            bool lmBaseFound=false;
            for (size_type r=0,h=i+1; r<getMaxLength(i+1); r++,h=rb(h)) {
                // leftmost base
                if (!lmBaseFound && baseIndex[h].first != -1) {
                    baseIndex[i].first=baseIndex[h].first;
                    lmBaseFound=true;
                }

                // rightmost base
                if (baseIndex[h].second != -1) {
                    baseIndex[i].second=baseIndex[h].second;
                }
            }

            // report pairing bases if P node
            if (c==ALPHA_BASEPAIR) {
                assert(baseIndex[i].first != -1);
                assert(baseIndex[i].second != -1);
                assert(baseIndex[i].first < baseIndex[i].second);

                pairs[baseIndex[i].first]=baseIndex[i].second;
                pairs[baseIndex[i].second]=baseIndex[i].first;
            }
        }
    }

    delete[] baseIndex;
}

/* ****************************************** */
/*             Public functions               */
/* ****************************************** */

void RNA_Alignment::getSequenceAlignments(std::string &s1, std::string &s2) const {
    s1="";
    s2="";

    // generate base strings
    for (size_type i=0; i<size(); i++) {
        if (lb_[i].a != ALPHA_BASEPAIR && lb_[i].b != ALPHA_BASEPAIR) {
            s1 += lb_[i].a;
            s2 += lb_[i].b;
        }
    }
}

void RNA_Alignment::getStructureAlignment(std::string &s, bool first) const {
    s="";

    std::map<unsigned int,unsigned int> pairs;
    makePairTable(pairs, first);

    // iterate through leaves nodes and use information of pairs
    for (size_type i=0; i<size(); i++) {
        RNA_Alphabet c;

        if (first)
            c=lb_[i].a;
        else
            c=lb_[i].b;

        if (isLeaf(i)) {
            if (c != ALPHA_GAP) {
                if (pairs.find(i) != pairs.end()) {	// is base paired ?
                    size_type j=pairs[i];
                    if (i<j)
                        s+='(';
                    else
                        s+=')';
                } else
                    s+='.';
            } else
                s+='-';
        }
    }
}

#ifdef HAVE_LIBG2
#ifdef HAVE_LIBRNA
void RNA_Alignment::squigglePlot(const std::string &filename_suffix, const RNAFuncs::SquigglePlotOptions &options) const {
    std::string str1,str2,seq1,seq2;
    std::string filename;

    getStructureAlignment(str1,true);
    getStructureAlignment(str2,false);
    getSequenceAlignments(seq1,seq2);

    filename = "x_" + filename_suffix;
    RNAFuncs::drawRNAAlignment(str1,str2,seq1,seq2,strname1_,strname2_,filename,true,options);
    filename = "y_" + filename_suffix;
    RNAFuncs::drawRNAAlignment(str2,str1,seq1,seq2,strname1_,strname2_,filename,false,options);
}
#endif
#endif

void RNA_Alignment::setStructureNames(const std::string &s1,const std::string &s2) {
    strname1_=s1;
    strname2_=s2;
}

#ifdef HAVE_LIBRNA
void RNA_Alignment::generateXML(std::ostream &s) const {
    std::string str1,str2,seq1,seq2;

    getStructureAlignment(str1,true);
    getStructureAlignment(str2,false);
    getSequenceAlignments(seq1,seq2);

    RNAFuncs::generateRNAAlignmentXML(str1,str2,seq1,seq2,strname1_,strname2_,s);
}
#endif


