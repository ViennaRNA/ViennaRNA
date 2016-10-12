#ifndef PARSER_H
#define PARSER_H

#include <string>
#include "rnafuncs.h"
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
#include <libxml++/libxml++.h>
#include <libxml/xmlschemas.h>
//#include "xsd.h"
#endif
#endif


std::istream * parseXMLFile(std::string filename, std::istream *inputStream, RNAFuncs::AddXmlInfos & xmlInfos);

void cutAfterChar(std::string &s,char c);

void structure2anchoredStructure(const std::string &sequence,
		    const std::string &structure, const std::string &shape,
				    std::string &anchoredViennaStructure);
void parseError(std::string input_type, int line_nr, std::string message); 

#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
extern "C" {
    bool validateXSD(std::string filename);
}
#endif
#endif

#endif  // PARSER_H 
