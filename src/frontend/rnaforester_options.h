#ifndef _RNAFORESTER_OPTIONS_H
#define _RNAFORESTER_OPTIONS_H

//#ifndef WIN32
//#include "config.h"
//#endif

//#ifndef HAVE_LIBRNA
//#undef HAVE_LIBG2
//#endif

#include "config.h"
#include "Arguments.h"

class Options {
private:
    struct OptionsInfo {
        std::string tag;
        std::string parameter;
        std::string filler;
        std::string description;
        bool hidden;
    };

public:
    enum RNAforesterOption {
			  // general options
        Help=0,
        SecretHelp,
        Version,
        ReadFromFile,
        SpaceTimeInfo,
        ShowOnlyScore,
        NoScale,
        Tables,
        Backtrace,
				// computation properties
				Topdown,
        CalculateDistance,
        RelativeScore,
        LocalSimilarity,
        LocalSubopts,
        SmallInLarge,
				Anchoring,
				Affine,
        Multiple,
        RIBOSUMScore,
        TreeEdit,
        GlobalAlignment,
				// progressive alignment
        ConsensusMinBaseProb,
        ConsensusMinPairProb,
        ClusterThreshold,
        ClusterJoinCutoff,
#ifdef HAVE_LIBRNA
        PredictProfile,
        PredictMinPairProb,
#endif
        SaveProfile,
        ProfileSearch,
				// score values
        BpRepScore,
        BpDelScore,
        BpDelOpenScore,
        BMatchScore,
        BRepScore,
        BDelScore,
        BDelOpenScore,
				// squiggle plot options
#ifdef HAVE_LIBG2
        MakeSquigglePlot,
        SquiggleHideBaseNumbers,
        SquiggleBaseNumberInterval,
        SquiggleGreyColors,
        SquiggleScaleFactor,
        SquiggleGenerateFIG,
#ifdef HAVE_LIBGD
        SquiggleGeneratePNG,
        SquiggleGenerateJPG,
#endif
#endif
				// other output options
        FastaOutput,
        MakeDotForInputTrees,
        OutputAlignmentDotFormat,
#ifdef HAVE_LIBXMLPLUSPLUS
#ifdef HAVE_LIBXML2
        GenerateXML,
        XmlOutputFile,
#endif
#endif
        NumberOfOptions
    };

    class IncompatibleException {
    private:
        std::string tag1_;
        std::string tag2_;

    public:
        IncompatibleException(std::string tag1,std::string tag2)
                : tag1_(tag1), tag2_(tag2) {};

        void showError();
    };

    class RequiresException {
    private:
        std::string tag1_;
        std::string tag2_;

    public:
        RequiresException(std::string tag1,std::string tag2)
                : tag1_(tag1), tag2_(tag2) {};

        void showError();
    };

    Options(int argc, const char **argv);
    ~Options();

    bool has(RNAforesterOption option) const;
    const char** getArgs() const;
    const unsigned int getNrOfOptions() const;

    template<class T>
    void get(RNAforesterOption option, T &var, T def) const {
        args_->get(options_[option].tag,var,def);
    }

    void help();
    void secretHelp();
    std::string generateFilename(RNAforesterOption option, const std::string &suffix, const std::string &defName, unsigned int count=0) const;

private:
    OptionsInfo *options_;
    Arguments *args_;
    const char **argv_;
    unsigned int nrArgs;

    inline void setOption(RNAforesterOption,std::string tag, std::string parameter, std::string filler, std::string description, bool hidden);
    void exclude(RNAforesterOption opt1, RNAforesterOption opt2);
    void requires(RNAforesterOption opt1, RNAforesterOption opt2);

};

#endif
