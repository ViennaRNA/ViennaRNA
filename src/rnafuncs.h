#ifndef _RNA_FUNCS_H
#define _RNA_FUNCS_H

#include <vector>
#include <string>
#include "rnaforester_options.h"
#include "rna_profile_alignment.h"

#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
#include <libxml++/libxml++.h>
#endif
#endif



class RNAFuncs {
public:
    struct SquigglePlotOptions {
        bool hideBaseNumbers;
        unsigned int baseNumInterval;
        bool greyColors;
        bool generatePNG;
        bool generateJPG;
        bool generateFIG;
        double scale;
    };

    struct AddXmlInfos {
        std::map<int,std::string> idmapping;
        std::map<int,std::string> comments;
        std::map<int,std::string> descriptions;
        std::map<int,std::string> names;
        std::map<int,std::string> synonyms;
        unsigned int xbasepos,ybasepos;
        bool xmlInput;
    };

    static bool isRNAString(const std::string &str);
    static bool isViennaString(const std::string &str, unsigned long &basePairCount, unsigned long &maxDepth);
		static bool isViennaString(const std::string &string, unsigned long &length, unsigned long &basePairCount, unsigned long &maxDepth);
		static bool isAnchoredViennaString(const std::string &string, unsigned long &length, unsigned long &basePairCount, unsigned long &maxDepth);

    static void drawRNAStructure(const std::string &seq, const std::string &structure, const std::string &filename_prefix, const std::string &structname, const std::vector<std::pair<unsigned int,unsigned int> > &regions, const SquigglePlotOptions &options);
    static void drawRNAAlignment(const std::string &structure, const std::string &altStructure,  const std::string &seq1, const std::string &seq2, const std::string &strname1, const std::string &strname2, const std::string &filename_prefix, const bool atX, const SquigglePlotOptions &options);
    static void generateRNAAlignmentXML(const std::string &structure, const std::string &altStructure, const std::string &seq1, const std::string &seq2, const std::string &strname1, const std::string &strname2, std::ostream &s);
    static void printAli(const std::string &name1, const std::string &id2, const std::string &seq1, const std::string &seq2, const std::string &str1, const std::string &str3);

    static unsigned int treeSize(const std::string &viennaStr, bool anchored);

#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
    static std::string getXSDURL();
    static void printMAliXML(std::vector<std::pair<double,RNAProfileAlignment*> > &resultList, const Options &options,double &minPairProb,AddXmlInfos &xmlInfos,const std::string &outputFile);
    static void printPAliXML(const std::string &id1, const std::string &name2, const std::string &seq1, const std::string &seq2, const std::string &str1, const std::string &str3, double &score, const Options &options,AddXmlInfos &xmlInfos,const std::string &outputFile);
    static void printMapping(std::map<int,std::string> &mapping);
#endif
#endif
    static std::string UpperCase(const std::string &str);
};

#endif




