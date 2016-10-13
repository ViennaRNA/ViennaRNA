#include <iostream>
#include <sstream>

//#ifndef WIN32
//#include "config.h"
//#endif

#include "config.h"
//#ifndef HAVE_LIBRNA
//#undef HAVE_LIBG2
//#endif

#include "rnaforester_options.h"

Options::Options(int argc, const char **argv)
        : argv_(argv) {
    nrArgs = argc;
    options_=new OptionsInfo[NumberOfOptions];

		// general options
    setOption(Help,                      "--help","","                    ","shows this help info",false);
    setOption(Version,                   "--version","","                 ","shows version information",false);
    setOption(OutputAlignmentDotFormat,  "-dot","=file","                 ","show alignment forest in dot format",true);
    setOption(Tables,                    "--tables","","                  ","shows dynamic programming tables",false);
    setOption(Backtrace,                 "--backtrace","","               ","shows backtrace call table cells",false);
		// computation properties
    setOption(Topdown,					         "-t","","                        ","calculate alignment top down instead of bottom up",false);
    setOption(CalculateDistance,         "-d","","                        ","calculate distance instead of similarity",false);
    setOption(RelativeScore,             "-r","","                        ","calculate relative score",false);
    setOption(LocalSimilarity,           "-l","","                        ","local similarity",false);
    setOption(LocalSubopts,		           "-so","=int","                   ","local suboptimal alignments within int%",false);
    setOption(SmallInLarge,              "-s","","                        ","small-in-large similarity",false);
    setOption(Anchoring,                 "--anchor","","                  ","use shape anchoring for speedup",false);
    setOption(Affine,                    "-a","","                        ","affine gap scoring",false);
    setOption(Multiple,                  "-m","","                        ","multiple alignment mode",false);
    setOption(ClusterThreshold,          "-mt","=double","                ","clustering threshold",false);
    setOption(ClusterJoinCutoff,         "-mc","=double","                ","clustering cutoff",false);
#ifdef HAVE_LIBRNA
    setOption(PredictProfile,            "-p","","                        ","predict structures from sequences",false);
    setOption(PredictMinPairProb,	       "-pmin","=double","              ","minimum basepair frequency for prediction",false);
#endif
    setOption(SaveProfile,	             "-sp","=file","                  ","save profile",false);
    setOption(ProfileSearch,	           "-ps","=file","                  ","profile search",false);

    setOption(TreeEdit,                  "-e","","                        ","use tree edit model to calculate distance and similarity",true);
    setOption(GlobalAlignment,           "-g","","                        ","calculate global alignment in its original version",true);
		// score properties
    setOption(BpRepScore,                "-pm","=int","                   ","basepair(bond) match score",false);
    setOption(BpDelOpenScore,            "-pdo","=int","                  ","basepair bond indel open score",false);
    setOption(BpDelScore,                "-pd","=int","                   ","basepair bond indel score",false);
    setOption(BMatchScore,               "-bm","=int","                   ","base match score",false);
    setOption(BRepScore,                 "-br","=int","                   ","base mismatch score",false);
    setOption(BDelOpenScore,             "-bdo","=int","                  ","base indel open score",false);
    setOption(BDelScore,                 "-bd","=int","                   ","base indel score",false);
    setOption(RIBOSUMScore,              "--RIBOSUM","","                 ","RIBOSUM85-60 scoring matrix (base-pair substitutions)",false);
    setOption(ConsensusMinBaseProb,	     "-cbmin","=double","             ","minimum base frequency for consensus structure",false);
    setOption(ConsensusMinPairProb,	     "-cmin","=double","              ","minimum basepair frequency for consensus structure",false);
		// squiggleplot options
#ifdef HAVE_LIBG2
    setOption(MakeSquigglePlot,          "-2d","","                       ","generate alignment 2D plots in postscript format",false);
    setOption(SquiggleHideBaseNumbers,   "--2d_hidebasenum","","          ","hide base numbers in 2D plot",false);
    setOption(SquiggleBaseNumberInterval,"--2d_basenuminterval","=n","    ","show every n-th base number",false);
    setOption(SquiggleGreyColors,        "--2d_grey","","                 ","use only grey colors in 2D plots",false);
    setOption(SquiggleScaleFactor,       "--2d_scale","=double","         ","scale factor for the 2d plots",false);
    setOption(SquiggleGenerateFIG,       "--2d_fig","","                  ","generate additional fig file of 2d plot",false);
#ifdef HAVE_LIBGD
    setOption(SquiggleGeneratePNG,       "--2d_png","","                  ","generate additional png file of 2d plot",false);
    setOption(SquiggleGenerateJPG,       "--2d_jpg","","                  ","generate additional jpg file of 2d plot",false);
#endif
#endif
    setOption(ReadFromFile,              "-f","=file","                   ","read input from file",false);
    //  setOption(SaveMultipleAliFile,       "-sm","=file","                  ","save multiple alignment as binary file",true);
    setOption(NoScale,                   "--noscale","","                 ","suppress output of scale",false);
    setOption(MakeDotForInputTrees,      "-idot","","                     ","make dot files for the input trees",true);
#ifdef HAVE_LIBRNA
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
    setOption(GenerateXML,               "--xml","","                     ","generate xml file in RNAStructAlignmentML format",false);
    setOption(XmlOutputFile,             "--xml_output","","                     ","name of xml output file",false);
#endif
#endif
#endif
    setOption(SecretHelp,                "--shelp","","                   ","shows this help info",true);
    setOption(ShowOnlyScore,             "--score","","                   ","compute only scores, no alignment",false);
    setOption(FastaOutput,               "--fasta","","                   ","generate fasta output of alignments",false);
    setOption(SpaceTimeInfo,             "--spacetime","","               ","space and time measurements",true);

    // set the arguments that can be seperated by spaces
    std::stringstream ss;
    for (int i=0; i<NumberOfOptions; i++) {
        if (!options_[i].parameter.empty())
            ss << options_[i].tag << "|";
    }

    Arguments::setArgumentsWithSpaces(ss.str());
    args_ = new Arguments(argc,argv);

    // check options for compatibility
    exclude(LocalSimilarity,SmallInLarge);
#ifdef HAVE_LIBG2
    exclude(SquiggleHideBaseNumbers,Multiple);
    exclude(SquiggleBaseNumberInterval,Multiple);
#endif
    exclude(CalculateDistance,LocalSimilarity);
    exclude(CalculateDistance,RIBOSUMScore);
    exclude(CalculateDistance,RelativeScore);
    exclude(CalculateDistance,Multiple);
    exclude(CalculateDistance,SmallInLarge);
    exclude(Multiple,RIBOSUMScore);
    exclude(RIBOSUMScore,BpRepScore);
    exclude(RIBOSUMScore,BMatchScore);
    exclude(RIBOSUMScore,BDelScore);
    exclude(LocalSimilarity,Topdown);

		requires(Anchoring, Topdown);
    requires(LocalSubopts,LocalSimilarity);
    requires(ClusterThreshold,Multiple);
    requires(ClusterJoinCutoff,Multiple);
#ifdef HAVE_LIBRNA
    requires(PredictProfile,Multiple);
    requires(PredictMinPairProb,PredictProfile);
#endif
}

Options::~Options() {
    delete[] options_;
    delete args_;
}

inline void Options::setOption(RNAforesterOption i,std::string tag, std::string parameter, std::string filler, std::string description, bool hidden) {
    options_[i].tag=tag;
    options_[i].parameter=parameter;
    options_[i].filler=filler;
    options_[i].description=description;
    options_[i].hidden=hidden;
}

bool Options::has(RNAforesterOption option) const {
    return args_->has(options_[option].tag);
}

const char** Options::getArgs() const {
    return argv_;
}

const unsigned int Options::getNrOfOptions() const {
    return nrArgs;
}

void Options::help() {
    std::cout << "Usage: " << argv_[0] << " [options]" << std::endl;

    for (int i=0; i<NumberOfOptions; i++) {
        if (!options_[i].hidden)
            std::cout << options_[i].tag << options_[i].parameter << options_[i].filler << options_[i].description << std::endl;
    }
}

void Options::secretHelp() {
    std::cout << "Usage: " << argv_[0] << " [options]" << std::endl;
    std::cout << "help including hidden parameters" << std::endl;

    for (int i=0; i<NumberOfOptions; i++) {
        std::cout << options_[i].tag << options_[i].parameter << options_[i].filler << options_[i].description;
        if (options_[i].hidden)
            std::cout << " *";
        std::cout << std::endl;
    }
}

std::string Options::generateFilename(RNAforesterOption option, const std::string &suffix, const std::string &defName, unsigned int count) const {
    std::string s;

    get(option,s,std::string(""));
    if (s=="") {
        if (has(ReadFromFile)) {
            std::ostringstream ss;

            get(ReadFromFile,s,defName);
            ss << s;
            if (count)
                ss << "_" << count;

            ss << suffix << '\0';

            s = ss.str();
        } 
				else {
            s = defName;
        }
    }
    return s;
}

void Options::exclude(RNAforesterOption opt1, RNAforesterOption opt2) {
    if (has(opt1) && has(opt2))
        throw IncompatibleException(options_[opt1].tag,options_[opt2].tag);
}

void Options::requires(RNAforesterOption opt1, RNAforesterOption opt2) {
    if (has(opt1) && !has(opt2))
        throw RequiresException(options_[opt1].tag,options_[opt2].tag);
}


void  Options::IncompatibleException::showError() {
    std::cerr << "The options " << tag1_ << " and " << tag2_ << " exclude each other." << std::endl;
}

void  Options::RequiresException::showError() {
    std::cerr << "The option " << tag1_ << " requires option " << tag2_ << "." <<  std::endl;
}
