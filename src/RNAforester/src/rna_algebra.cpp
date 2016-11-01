#include "rna_algebra.h"

Score::Score(Options &options) 
	: isAffine_(options.has(Options::Affine)),
		isLocal_(options.has(Options::LocalSimilarity)),
		isDistance_(options.has(Options::CalculateDistance)),
		isRIBOSUM_(options.has(Options::RIBOSUMScore)),
    bp_rep_score_(0),
    bp_indel_score_(0),
    bp_indel_open_score_(0),
    b_match_score_(0),
    b_rep_score_(0),
    b_indel_score_(0),
    b_indel_open_score_(0) {



    // distance or similarity ?
    if (isDistance_) {
        bp_rep_score_  = 0;
        bp_indel_score_ = 3;
        b_match_score_ = 0;
        b_rep_score_   = 1;
        b_indel_score_ = 2;
				if (isAffine_) {
					bp_indel_open_score_ = 3; // TODO compl. arbitrary!
					b_indel_open_score_ = 2; // TODO compl. arbitrary!;
				}
    } 
		else {
        if (isRIBOSUM_) {
            bp_indel_score_ =-100;
            b_indel_score_  =-200;
        } 
				else {
            bp_rep_score_  = 10;
            bp_indel_score_ =-5;
            b_match_score_ = 1;
            b_rep_score_   = 0;
            b_indel_score_ =-10;
						if (isAffine_) {
							bp_indel_open_score_ = -5; // TODO compl. arbitrary!
							b_indel_open_score_ = -10; // TODO compl. arbitrary!;
						}
        }
    }

    // read scores
    options.get(Options::BpRepScore,bp_rep_score_,bp_rep_score_);
    options.get(Options::BpDelScore,bp_indel_score_,bp_indel_score_);
    options.get(Options::BMatchScore,b_match_score_,b_match_score_);
    options.get(Options::BRepScore,b_rep_score_,b_rep_score_);
    options.get(Options::BDelScore,b_indel_score_,b_indel_score_);
		if (isAffine_) {
			options.get(Options::BpDelOpenScore,bp_indel_open_score_,bp_indel_open_score_);
			options.get(Options::BDelOpenScore,b_indel_open_score_,b_indel_open_score_);
		}
}

Score::Score(const Score &s) {
    // copy scores
		isAffine_ = s.isAffine_;
		isLocal_ = s.isLocal_;
    isDistance_ = s.isDistance_;
    isRIBOSUM_ = s.isRIBOSUM_;
    bp_rep_score_ = s.bp_rep_score_;
    bp_indel_score_ = s.bp_indel_score_;
    bp_indel_open_score_ = s.bp_indel_open_score_;
    b_match_score_ = s.b_match_score_;
    b_rep_score_ = s.b_rep_score_;
    b_indel_score_ = s.b_indel_score_;
    b_indel_open_score_ = s.b_indel_open_score_;
}

std::ostream& operator<< (std::ostream &out, const Score &score) {

    out << "*** Scoring parameters ***" << std::endl << std::endl;

    out << "Scoring type: ";
    if (score.isAffine_)
        out << "affine ";
    if (score.isDistance_)
        out << "distance" << std::endl;
    else {
        if (score.isLocal_)
            out << "local ";
				else
            out << "global ";
        out << "similarity" << std::endl;
    }

    std::cout << "Scoring parameters:" << std::endl;
    if (score.isRIBOSUM_) {
        out << "RIBOSUM85-60 Scoring matrix" << std::endl;
        out << "pair indel:  " << score.bp_indel_score_ << std::endl;
        out << "base indel:  " << score.b_indel_score_ << std::endl << std::endl;
    } 
		else {
        out << "pair match:       " << score.bp_rep_score_ << std::endl;
				if (score.isAffine_)
					out << "pair indel open:  " << score.bp_indel_open_score_ << std::endl;
        out << "pair indel:       " << score.bp_indel_score_ << std::endl;
        out << "base match:       " << score.b_match_score_ << std::endl;
        out << "base replacement: " << score.b_rep_score_ << std::endl;
        out << "base indel:       " << score.b_indel_score_ << std::endl;
				if (score.isAffine_)
					out << "base indel open:	" << score.b_indel_open_score_ << std::endl;
        out << std::endl;
    }


	return out;
}
