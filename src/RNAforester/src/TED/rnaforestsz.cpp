#include "rnaforestsz.h"
#include "rnafuncs.h"

#include "forestsz.t.cpp"

RNAForestSZ::RNAForestSZ(const std::string &baseStr, const std::string &viennaStr, const std::string &name)
        : ForestSZ<RNA_Alphabet>(RNAFuncs::treeSize(viennaStr,false)),
        baseStr_(baseStr), viennaStr_(viennaStr) {
    unsigned int pos,node;
    unsigned long basePairCount,maxDepth;

    assert(RNAFuncs::isRNAString(baseStr));
    RNAFuncs::isViennaString(viennaStr,basePairCount,maxDepth);

    // check if base string and vienna string have the same length
//  if(baseStr.length() > 0)
//    if(baseStr.length() != viennaStr_.length())
//      throw RNAForestExceptionInput(RNAForestExceptionInput::Error_BaseStringAndViennaStringIncompatible);

    if (name.empty())
        name_="unknown";
    else
        name_=name;

    pos=0;
    node=0;
    buildForest(pos,node);
    calcKeyroots();

    // debug
    /*
    int i;
    std::cout << "lb_:" << std::endl;
    for(i=0;i<size_;i++)
        std::cout << i << ": " << RNA_Alpha2alpha(lb_[i]) << std::endl;

    std::cout << "lml_:" << std::endl;
    for(i=0;i<size_;i++)
        std::cout << i << ": " << lml_[i] << std::endl; */
}

void RNAForestSZ::buildForest(unsigned int &pos, unsigned int &node) {
    unsigned int node2;

    switch (viennaStr_[pos]) {
    case '.':
        // set label
        if (baseStr_.length())
            lb_[node]=alpha2RNA_Alpha(baseStr_[pos]);
        else
            lb_[node]=ALPHA_BASE;

        lml_[node]=node;

        break;
    case '(':
        // left paring base
        if (baseStr_.length())
            lb_[node]=alpha2RNA_Alpha(baseStr_[pos]);
        else
            lb_[node]=ALPHA_BASE;

        lml_[node]=node;

        // build subforest right to (
        node2=node;
        node++;
        pos++;
        buildForest(pos,node);


        // right bracket
        assert(viennaStr_[pos]==')');

        // right pairing base
        if (baseStr_.length())
            lb_[node]=alpha2RNA_Alpha(baseStr_[pos]);
        else
            lb_[node]=ALPHA_BASE;

        lml_[node]=node;


        node++;
        lb_[node]=ALPHA_BASEPAIR;
        lml_[node]=node2;
        break;

    case ')':
        return;
        break;
    }

    pos++;
    node++;
    if (pos<viennaStr_.length())
        buildForest(pos,node);
}
