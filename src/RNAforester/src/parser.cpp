#include "parser.h"
#include "anchors.h"


std::istream * parseXMLFile(std::string filename, std::istream *inputStream, RNAFuncs::AddXmlInfos & xmlInfos) {
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
    // create Dom parser
    xmlpp::DomParser domParser;
    domParser.parse_file(filename);

    xmlpp::Document* document = domParser.get_document();
    xmlpp::Element* rootElement = document->get_root_node();

    int seqID = 0;
    std::string foresterFormat = "";
    // temporary maps for xmlinfo
    std::map<int,std::string> comments, idMapping, descriptions, names, synonyms;

    // for each rnastructure element
    xmlpp::Node::NodeList::iterator it1, it2, it3;

    xmlpp::Node::NodeList structureNodes = rootElement->get_children();
    for (it1=structureNodes.begin(); it1!=structureNodes.end(); it1++) {
        if ((*it1)->get_name() == "rnastructure") {
            seqID++;
            std::stringstream ss;
            ss << ">" << seqID << "\n";
            foresterFormat.append(ss.str());
            ss.str("");
        }

        xmlpp::Node::NodeList sequenceNodes = (*it1)->get_children();
        for (it2 = sequenceNodes.begin(); it2!= sequenceNodes.end(); it2++) {
            if ((*it2)->get_name()=="sequence") {
                xmlpp::Element* elemSeq = (xmlpp::Element*)(*it2);

                // map internal id to seqID
                std::cout << "seqID: " << seqID << std::endl;
                std::cout << "idMapping[seqID]: " << elemSeq->get_attribute("seqID")->get_value() << std::endl;
                idMapping[seqID] = elemSeq->get_attribute("seqID")->get_value();

                xmlpp::Node::NodeList nameNodes = elemSeq->get_children();
                for (it3=nameNodes.begin(); it3!=nameNodes.end(); it3++) {
                    if ((*it3)->get_name()=="name" && ((xmlpp::Element*)(*it3))->get_child_text()) {
                        // map the name to seqID
                        names[seqID] = ((xmlpp::Element*)(*it3))->get_child_text()->get_content();
                    }
                    if ((*it3)->get_name()=="synonyms" && ((xmlpp::Element*)(*it3))->get_child_text()) {
                        // map synonyms to seqID
                        synonyms[seqID] = ((xmlpp::Element*)(*it3))->get_child_text()->get_content();
                    }
                    if ((*it3)->get_name()=="description" && ((xmlpp::Element*)(*it3))->get_child_text()) {
                        // map the description to seqID
                        descriptions[seqID] = ((xmlpp::Element*)(*it3))->get_child_text()->get_content();
                    }
                    if ((*it3)->get_name()=="freeSequence" || (*it3)->get_name()=="nucleicAcidSequence") {
                        xmlpp::Element* elemFreeSequence = (xmlpp::Element*)(*it3);
                        foresterFormat.append(elemFreeSequence->get_child_text()->get_content());
                        foresterFormat.append("\n");
                    }
                } // end for it3
            }

            // get comment if available
            if ((*it2)->get_name()=="comment" && ((xmlpp::Element*)(*it2))->get_child_text()) {
                comments[seqID] = ((xmlpp::Element*)(*it2))->get_child_text()->get_content();
            }

            // get structure
            if ((*it2)->get_name()=="structure" && ((xmlpp::Element*)(*it2))->get_child_text()) {
                xmlpp::Element* elemStr = (xmlpp::Element*)(*it2);
                foresterFormat.append(elemStr->get_child_text()->get_content());
                foresterFormat.append("\n");
            }
        } // end for it2

    } // end for it1

    // struct with additional xml infos
    xmlInfos.idmapping = idMapping;
    xmlInfos.comments = comments;
    xmlInfos.descriptions = descriptions;
    xmlInfos.names = names;
    xmlInfos.synonyms = synonyms;
    xmlInfos.xmlInput = true;

    //std::istringstream * inputString = new std::istringstream(foresterFormat);
    inputStream = new std::istringstream(foresterFormat);
    return inputStream;
#endif
#endif
}


void cutAfterChar(std::string &s,char c) {
    std::string::size_type pos=s.find(c);
    if (pos!=std::string::npos)
        s.erase(pos);
}


// transforms struct+shape into intermediate representation,
// that is a vienna struct with anchor symbols before anchored basepairs

void structure2anchoredStructure(const std::string &sequence,
    const std::string &structure, const std::string &shape,
    std::string &anchoredViennaStructure) {

  std::list<int> nodesWithAnchors;
	// call bison parser in anchors directory for parsing
  parseStructure(structure, &nodesWithAnchors);

  // check whether we found the correct anchors as in the given shape by rnacast
  unsigned int openBr = 0;
  unsigned int closeBr = 0;
  for (unsigned int i = 0; i < shape.size(); i++) {
    if (shape[i] == '[')
      openBr++;
    else if (shape[i] == ']')
      closeBr++;
  }
  if (openBr != closeBr)
    parseError("vienna (anchoring)", 0, "shape has different number of opening and    closing brackets");
  if (openBr != nodesWithAnchors.size())
    parseError("vienna (anchoring)", 0, "given shape brackets imply different anchors than a new folding into this shape");
  std::cout << "Tree with " << nodesWithAnchors.size() << " anchors." << std::endl;

  // add 'a' symbols before anchored basepairs in structure
  std::ostringstream oss;
  for (unsigned int i = 0; i < structure.size(); i++) {
    if ((nodesWithAnchors.size() != 0)
        && (nodesWithAnchors.front() == (int)i)) { // anchor intended at this position
      oss << "a" << structure[i];
      nodesWithAnchors.pop_front();
    }
    else {
      oss << structure[i];
    }
  }
  anchoredViennaStructure = oss.str();
  return;
}

void parseError(std::string input_type, int line_nr, std::string message) {
  std::cerr << "Error parsing input of type " << input_type;
  std::cerr << " in line " << line_nr << ":" << std::endl;
  std::cerr << message << "." << std::endl;
  exit(EXIT_FAILURE);
}

#ifdef HAVE_LIBXML2
#ifdef HAVE_LIBXMLPLUSPLUS
extern "C" {
    bool validateXSD(std::string filename) {
        xmlSchemaParserCtxtPtr ctxt;
        xmlSchemaValidCtxtPtr valCtxt;
        xmlSchemaPtr schema;
        int val;
        // TODO den gab es nicht!
        const char* XSD_STRING = "";

        ctxt = xmlSchemaNewMemParserCtxt(XSD_STRING,sizeof(XSD_STRING));

        schema = xmlSchemaParse(ctxt);

        valCtxt = xmlSchemaNewValidCtxt(schema);

        val = xmlSchemaValidateFile(valCtxt,filename.c_str(),0);

        if (val==0) {
            return true;
        } else {
            return false;
        }
    }
}
#endif
#endif

