#ifndef _RNA_FOREST_H_
#define _RNA_FOREST_H_

#include <iostream>

#include "algebra.h"
#include "forest.h"
#include "rna_alphabet.h"
#include "rnafuncs.h"





/** RNAForest encapsulates the forest representation of an RNA */
class RNAForest
            : public Forest<RNA_Alphabet> {
private:
    typedef Forest<RNA_Alphabet> ForestRNAType;

    std::string name_;
    std::string baseStr_;
    std::string viennaStr_;

    void showLabel(std::ostream &s,RNA_Alphabet a) const;	// virtual function	of Forest
    void makeLabel(RNA_Alphabet &a,char c);

    /** Build	forest from	structure, sequence	pair. */
    void buildForest(const std::string	&baseStr, const	std::string &viennaStr);
    void buildForestAnchored(const std::string	&baseStr, const	std::string &viennaStr);

public:
    /**	Build forest from structure	sequence pair. This allows direct use of Vienna RNA strings.
    *	e.g.: str="(..(...))", seq="accguuucu"
    */
    RNAForest(const std::string &baseStr, const std::string	&viennaStr, const std::string &name, bool anchored);
		virtual ~RNAForest();
    const	std::string&	getName() const	{
        return	name_;
    };
    const	std::string&	getBaseStr() const {
        return baseStr_;
    };
    const	std::string&	getViennaStr() const {
        return viennaStr_;
    };

    void plot2d(const std::string &filename_prefix, const std::vector<std::pair<unsigned int,unsigned int> > &regions, const RNAFuncs::SquigglePlotOptions &sqOptions) const;
};

#endif


