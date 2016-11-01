#ifndef _RNA_FORESTSZ_H_
#define _RNA_FORESTSZ_H_

#include <iostream>

#include "forestsz.h"
#include "rna_alphabet.h"




// RNAForestSZ

class RNAForestSZ : public ForestSZ<RNA_Alphabet> {
private:
    std::string name_;
    std::string baseStr_;
    std::string viennaStr_;

    void buildForest(unsigned int &pos, unsigned int &node);
public:
    RNAForestSZ(const std::string &baseStr, const std::string &viennaStr, const std::string &name);

    const std::string& getName() const {
        return name_;
    };
    const std::string& getBaseStr() const {
        return baseStr_;
    };
    const std::string& getViennaStr() const {
        return viennaStr_;
    };

};

#endif


