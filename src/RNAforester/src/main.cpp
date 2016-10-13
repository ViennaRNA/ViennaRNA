//#ifndef WIN32
//#include "config.h"
//#endif

#include "config.h"
#include "rnaforester_options.h"
#include "parser.h"
#include "aligner.h"

static const std::string AUTHOR = "Stefanie Schirmer, Matthias Hoechsmann ";
static const std::string RNAFORESTER_VERSION = "2.0";
static const std::string PROMPT = "Input string (upper or lower case); & to end for multiple alignments, @ to quit\n";
static const std::string SCALE = "\e[22;36m....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8\e[0;0m\n";


static void handleSimpleOptions(Options & options, const char * progName) {
    if (options.has(Options::Help)) {
        options.help();
        exit(EXIT_SUCCESS);
    }
    if (options.has(Options::SecretHelp)) {
        options.secretHelp();
        exit(EXIT_SUCCESS);
    }
    if (options.has(Options::Version)) {
        std::cout << progName << ", version " << RNAFORESTER_VERSION << std::endl;
        std::cout << "Copyright " << AUTHOR << " 2001-2011," << std::endl;
        std::cout << "sschirme@techfak.uni-bielefeld.de" << std::endl;
        exit(EXIT_SUCCESS);
    }
}

static std::istream * readFromFile(Options &options, std::istream *inputStream, RNAFuncs::AddXmlInfos& xmlInfos) {
    std::string filename;
    std::string line;
    options.get(Options::ReadFromFile,filename,std::string(""));
    // NOTE: dont forget to delete in main loop
    std::ifstream * inputFile = new std::ifstream(filename.c_str());
    if (inputFile->fail()) {
        std::cerr << "cannot open file: \"" << filename << "\"" << std::endl;
        exit(EXIT_FAILURE);
    }
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
    // extracts chars of first line until '\n' to check for xml file
    getline(*inputFile,line);

    // input file is an xml file
    if (line.find("<?xml",0)==0) {
        return parseXMLFile(filename, inputStream, xmlInfos);
    } 
		else {
        xmlInfos.xmlInput = false;
        // no xml file - write the first line from line back to stream
        for (int i = line.size(); i>=0; i--) {
            inputFile->putback(line[i]);
        }
#endif
#endif
        return inputFile;
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
    }
#endif
#endif
}

// TODO durch die exits wird der speicher nicht mehr aufgeraeumt; unschoene dopplung bei anchoring da es auch in options steckt
static void align(std::vector<RNAProfileAlignment*> & inputListMult, std::vector<RNAForest*> & inputListPW,
		std::vector<RNAForestSZ*> & inputListSZ, bool multipleAlign, bool anchoring, Score & score, Options & options, RNAFuncs::AddXmlInfos & xmlInfos) {
    // ***** multiple alignment
    if ((options.has(Options::Multiple) && multipleAlign) || (options.has(Options::ProfileSearch) && inputListMult.size()==2)) {
        alignMultiple(inputListMult,score,options,anchoring,xmlInfos);

        if (options.has(Options::ProfileSearch)) {
            std::string filename;

            options.get(Options::ProfileSearch,filename,std::string(""));
            //RNAProfileAlignment *rnaProfileAli = new RNAProfileAlignment(filename,options.has(Options::Anchoring));
            RNAProfileAlignment *rnaProfileAli = new RNAProfileAlignment(filename,anchoring);
            inputListMult.push_back(rnaProfileAli);
						return;
        } 
				else
					return;
    }

    // ***** pairwise alignment
    if (inputListPW.size()==2) {

        if (options.has(Options::GlobalAlignment)) {
            alignPairwiseSimple(inputListPW,score,options,anchoring);
        } 
				else {
						clock_t start, finish;
						start = clock();
            alignPairwise(inputListPW,score,options,anchoring,xmlInfos);
						finish = clock();
						//std::cerr << "cpucyclesPerCodeblock " << ( (finish - start) ) << std::endl;
				}
				return;
    }

    if (inputListSZ.size()==2) {
        editPairwise(inputListSZ,score,options,anchoring);
				return;
    }
}

static void matchFields(std::string line, int lineNr, std::string *  energy, std::string *structure, std::string *shape, std::string *rank) {
		// remove leading whitespace
    line = line.substr(line.find_first_not_of("\t "),line.size());

    // match multiple fields in this line
    if (line.find_first_not_of("-.0123456789") == std::string::npos) {
      parseError("rnacast", lineNr, "No energy found");
    }
    *energy = line.substr(0, line.find_first_not_of("-.0123456789"));
    line = line.substr(line.find_first_not_of("-.0123456789"), line.size());
    line = line.substr(line.find_first_not_of(" "), line.size());

    if (line.find_first_not_of("(.)") == std::string::npos) {
      parseError("rnacast", lineNr, "No structure found");
    }
    *structure = line.substr(0, line.find_first_not_of("(.)"));
    line = line.substr(line.find_first_not_of("(.) "), line.size());

    if (line.find_first_not_of("[_]") == std::string::npos) {
      parseError("rnacast", lineNr, "No shape found");
    }
    *shape = line.substr(0, line.find_first_not_of("[_]"));
    line = line.substr(line.find_first_not_of("[_] "), line.size());

    if (line.find_first_of("-0123456789") == std::string::npos) {
      parseError("rnacast", lineNr, "No rank found - did you create without -o f (obsolete)?");
    }
    *rank = line.substr(line.find_first_of("-0123456789"), line.find_last_of("-0123456789"));
    assert(!structure->empty());

}

// parse RNAcast Shape Block
static void parseShapeBlock(Options &options, std::istream *inputStream, int& lineNr, Score &score, std::vector<RNAForest*> & inputListPW,std::vector<RNAProfileAlignment*>& inputListMult, RNAFuncs::AddXmlInfos& xmlInfos, int suboptPercent, double minPairProb) {
  std::string line = "";
  std::string sequence = "", name = "";
  unsigned long size = 0, basePairCount = 0, maxDepth = 0;
  int structure_count = 1;

  // 1. read a line - first line with the 1) is already processed and gone
  for (;;) {
    std::string structure = "", energy = "", shape = "", rank = "";
    getline(*inputStream, line);
		lineNr++;

		if (line.find_first_not_of("\t ") != std::string::npos )
			line = line.substr(line.find_first_not_of("\t "),line.size());

    // delete '\r' at line end from non unix files
    if (line[line.size() - 1] == '\r')
      line.erase(line.size() - 1);

    // 3. check for name of structure
    if (not line.empty() && line[0] == '>') {
      //nameStr = line;
      // cut >
      std::string fromGT = line.substr(1, line.size());
      // 4. cut blanks
      if (!fromGT.empty())
        name = fromGT.substr(fromGT.find_first_not_of(' '), fromGT.size());
      // set number if no name is given
      if (name.empty()) {
        std::ostringstream oss;
        oss << "Structure " << structure_count;
        name = oss.str();
      }
      assert(!name.empty());
      continue;
    }

    // 6. check for base string
    if (not line.empty() && RNAFuncs::isRNAString(line)) {
      sequence = line;
      // 7. cut it etc
      // convert to small letters
      std::transform(sequence.begin(), sequence.end(), sequence.begin(), ::tolower);
      // t -> u
      std::replace(sequence.begin(), sequence.end(), 't', 'u');
      assert(!sequence.empty());
      continue;
    }

    // 7. match line with structure and infos
		if (not line.empty() && line.rfind("R =") != std::string::npos) {
			matchFields(line, lineNr, &energy, &structure, &shape, &rank);
			if (structure.empty() || energy.empty() || shape.empty() || rank.empty()) {
				std::cout << "Could not parse line " << lineNr << ":"  << std::endl;
				std::cout << line << std::endl;
				std::cout << "Did you create the input without -o f (which is obsolete)?" << std::endl;
				exit(EXIT_FAILURE);
			}
		}

    // 8. no name given (predict prof) - use a counter s.o.
    // 9. check for vienna string   // check for vienna string
    if (not line.empty() && RNAFuncs::isViennaString(structure, size, basePairCount, maxDepth)) {
      assert(!structure.empty());
    }// 10. no name, no struc, no seq: error
    // unknown string
    else if (not line.empty()) {
      parseError("rnacast", lineNr, "The input sequence is neither a name nor an RNA//DNA string   nor a structure in dot-bracket format");
    }
    if (not line.empty() && sequence.length() != structure.length()) {
      parseError("rnacast", lineNr, "Sequence and structure differ in length");
    }
    // 11. add structure to input list - TODO be careful when we have no name(multiple)
    // 12. add if not multiple

    if (not line.empty() && options.has(Options::Anchoring)) {
      std::string anchoredViennaStructure = "";
      structure2anchoredStructure(sequence, structure, shape, anchoredViennaStructure);
      structure = anchoredViennaStructure;
			//std::cout << "anch. struct " << anchoredViennaStructure << std::endl;
    }

		bool anchoring = (options.has(Options::Anchoring));
		std::vector<RNAForestSZ*> inputListSZ;

		if (options.has(Options::Multiple)) {

			if (line.empty() || inputStream->eof() || line[0] == '@' || line[0] == '&') {
				//align
				bool multipleAlign = true;
				align(inputListMult,inputListPW,inputListSZ,multipleAlign,anchoring,score,options,xmlInfos);
				return;
			}
			// anchored or unanchored tree is built automatically via constructor
			RNAProfileAlignment *rnaProfileAli = new RNAProfileAlignment(sequence, structure, name, anchoring);
			inputListMult.push_back(rnaProfileAli);

		}
		else {
			// anchored or unanchored tree is built automatically via constructor
			RNAForest * rnaForest = new RNAForest(sequence, structure, name, anchoring);
			inputListPW.push_back(rnaForest);

			if (inputListPW.size() == 2) {

				// read until end of shape block
				while (not line.empty()) {
					getline(*inputStream, line);
					lineNr++;
				}
				// align
				bool multipleAlign = false;
				align(inputListMult,inputListPW,inputListSZ,multipleAlign,anchoring,score,options,xmlInfos);
				if (inputStream->eof() || line[0] == '@' || line[0] == '&')
					exit(EXIT_SUCCESS);
				return;
			}
		}

    name = "";
    sequence = "";
    // increase struct count
    structure_count++;
  }

}

// TODO suboptPercent macht gar nix, minPairProb aber schon
static void parseInputStream(Options &options, std::istream *inputStream, Score &score, std::vector<RNAProfileAlignment*>& inputListMult, RNAFuncs::AddXmlInfos& xmlInfos, int suboptPercent, double minPairProb) {
    std::vector<RNAForest*> inputListPW;
    std::vector<RNAForestSZ*> inputListSZ;
    std::string baseStr, viennaStr, nameStr;
    unsigned long basePairCount, maxDepth;
    unsigned int forest_count = 1;
    bool multipleAlign = false;
    bool showScale = false;
		bool anchored = options.has(Options::Anchoring);

		for (int lineNr = 0; ; lineNr++) {
        std::string line;
        getline(*inputStream,line);

        if (inputStream->eof()) {
            //if (options.has(Options::Multiple) && !options.has(Options::ProfileSearch))
            //    line="&";
            //else
                exit(EXIT_SUCCESS);
        }

        if (line.empty())
            continue;

        // quit if character is @
        if (line[0]=='@')
            break;

        // delete '\r' at line end from non unix files
        if (line[line.size()-1]=='\r')
            line.erase(line.size()-1);

        // check for name of structure
        if (line[0]=='>') {
            nameStr=&line[1];
            continue;
        }

        // cut after blank, but not for rnacast input
				if (line.rfind("R =") == std::string::npos && line.rfind("Shape") == std::string::npos )
					cutAfterChar(line,' ');

        // check for aligning multiple structures
        // if input is read from file the eof has the same meaning as &
        if ( line[0]=='&')
            multipleAlign = true;
        else {
            // check for base string

            if (RNAFuncs::isRNAString(line)) {
                baseStr=line;
                // convert to small letters
                transform(baseStr.begin(),baseStr.end(),baseStr.begin(),::tolower);
                // t -> u
                replace(baseStr.begin(),baseStr.end(),'t','u');
                // delete '.'  and '-' from alignment files
                remove(baseStr.begin(),baseStr.end(),'.');
                remove(baseStr.begin(),baseStr.end(),'-');

#ifdef HAVE_LIBRNA
                if (options.has(Options::PredictProfile)) {
                    std::ostringstream ss;
                    std::string constraint;

                    // if there is no name given (> ...) use a counter
                    if (nameStr=="")
                        ss << "> " << forest_count;
                    else
                        ss << nameStr;

                    std::cout << "Predicting structure profile for sequence: " << ss.str() << std::endl;
										// HERE NEW PROF ALI FOREST
                    RNAProfileAlignment *rnaProfileAli = new RNAProfileAlignment(baseStr,ss.str(),constraint,minPairProb);
                    inputListMult.push_back(rnaProfileAli);
                    makeDotFileInp(*rnaProfileAli,options,forest_count);
                    forest_count++;
                }
#endif

                //				  continue;
            } else {

                // check for vienna string
                if ((line.rfind("R =") == std::string::npos) && 
										(line.rfind("Shape") == std::string::npos) && 
										RNAFuncs::isViennaString(line,basePairCount,maxDepth)) {
									if (options.has(Options::Anchoring)) {
										std::cout << "Anchoring not possible without RNAcast input.\n";
										exit(EXIT_FAILURE);
									}

#ifdef HAVE_LIBRNA
                    // skip structure lines if structures are predicted
                    if (options.has(Options::PredictProfile)) {
                        std::cout << "ignoring structure: " << line << std::endl;
                        continue;
                    }
#endif
                    viennaStr = line;
                } 
								else if ((line.rfind("R =") != std::string::npos)
										|| (line.rfind("Shape") != std::string::npos)) {

									std::cout << "Parsing shape block " << line << std::endl;
									parseShapeBlock(options, inputStream, lineNr, score, inputListPW, inputListMult, xmlInfos, suboptPercent, minPairProb);
									if (inputStream->eof() || line[0] == '@' || line[0] == '&')
										exit(EXIT_SUCCESS);
									continue;
								}
								else {
                    std::cerr << "The input sequence is neither an RNA/DNA string nor in vienna format." << std::endl;
                    std::cerr << "line: " << line << std::endl;
                    showScale = true;
                    exit(EXIT_FAILURE);
                }


                // add structure to input vector
                if (options.has(Options::Multiple)) {
                    std::ostringstream ss;

                    // if there is no name given (> ...) use a counter
                    if (nameStr=="")
                        ss << "> " << forest_count;
                    else
                        ss << nameStr;

										// HERE NEW PROF ALI FOREST
                    RNAProfileAlignment *rnaProfileAli = new RNAProfileAlignment(baseStr,viennaStr,ss.str(),anchored);
                    makeDotFileInp(*rnaProfileAli,options,forest_count);
                    inputListMult.push_back(rnaProfileAli);
                } 
								else if (options.has(Options::TreeEdit)) {
                    RNAForestSZ *rnaForestSZ = new RNAForestSZ(baseStr,viennaStr,nameStr);
                    inputListSZ.push_back(rnaForestSZ);
                } 
								else {
                    RNAForest *rnaForest = new RNAForest(baseStr,viennaStr,nameStr,anchored);
										//rnaForest->printMembers();
                    nameStr = "";
                    makeDotFileInp(*rnaForest,options,forest_count);
                    inputListPW.push_back(rnaForest);
                }

                forest_count++;
                showScale = true;
            }
        }

				//std::cout << "align as usual." << std::endl;
				bool anchoring = false;
				align(inputListMult,inputListPW,inputListSZ,multipleAlign,anchoring,score,options,xmlInfos);
        multipleAlign = false;
        forest_count = 1;
    }
    // free dynamic allocated memory
    std::vector<RNAForest*>::const_iterator it;
    for (it = inputListPW.begin(); it!=inputListPW.end(); it++) {
        delete *it;
    }

}



int main(int argc, const char **argv) {

    try {
        // handle options
        Options options(argc,argv);
        handleSimpleOptions(options, argv[0]);

        // read and print score values
        Score score(options);
        if (!options.has(Options::ShowOnlyScore))
            std::cout << score;

        // set up input stream
        std::istream *inputStream = &std::cin;
        RNAFuncs::AddXmlInfos xmlInfos;
        if (options.has(Options::ReadFromFile)) {
            inputStream = readFromFile(options, inputStream, xmlInfos);
        }

        // check params for subopt and minpairprob
        int suboptPercent = 100;
        if (options.has(Options::LocalSubopts)) {
            options.get(Options::LocalSubopts,suboptPercent,100);
            if (suboptPercent<0 || suboptPercent>100) {
                std::cerr << "error: value for parameter --subopts must be in range from 0 to 100" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
        double minPairProb = 0.25;
#ifdef HAVE_LIBRNA
        options.get(Options::PredictMinPairProb,minPairProb,0.25);
        if (options.has(Options::PredictProfile))
            std::cout << "Minimum required basepair probability (-pmin): " << minPairProb << std::endl;
#endif
        if (!options.has(Options::ShowOnlyScore))
            if (options.has(Options::LocalSimilarity) && options.has(Options::LocalSubopts))
                std::cout << "calculate suboptimals within " << suboptPercent << "% of global optimum" << std::endl << std::endl;

        // profile search
        std::vector<RNAProfileAlignment*> inputListMult;
        if (options.has(Options::ProfileSearch)) {
            std::string filename;
            options.get(Options::ProfileSearch,filename,std::string(""));
            if (filename=="") {
                std::cerr << "no profile filename" << std::endl;
                exit(EXIT_FAILURE);
            }

            RNAProfileAlignment *rnaProfileAli = new RNAProfileAlignment(filename,options.has(Options::Anchoring));
            inputListMult.push_back(rnaProfileAli);
        }

        if (!options.has(Options::NoScale) && !options.has(Options::ReadFromFile))
            std::cout << std::endl << PROMPT << SCALE;

        // parse the input stream
        parseInputStream(options, inputStream, score, inputListMult, xmlInfos, suboptPercent, minPairProb);
        if (options.has(Options::ReadFromFile)) 
					delete inputStream; 
        return (0);
    }
    catch (Options::IncompatibleException e) {
        e.showError();
        return(EXIT_FAILURE);
    } 
		catch (Options::RequiresException e) {
        e.showError();
        return(EXIT_FAILURE);
    }
		catch (std::bad_alloc e) {
				std::cout << "Alignment computation aborted, out of memory." << std::endl;
        return(EXIT_FAILURE);
		}
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
    catch (xmlpp::validity_error ve) {
        std::cout << ve.what() << std::endl;
        //return(EXIT_FAILURE);
    } 
		catch (xmlpp::parse_error pe) {
        std::cout << pe.what() << std::endl;
        //return(EXIT_FAILURE);
    } 
		catch (xmlpp::exception e) {
        std::cout << e.what() << std::endl;
        //return(EXIT_FAILURE);
    }
#endif
#endif

}

