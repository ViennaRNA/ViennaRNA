#ifndef _RNA_ALGEBRA_H_
#define _RNA_ALGEBRA_H_

#include <assert.h>
#include <algorithm>
#include <climits>

#include "algebra.h"
#include "debug.h"
#include "misc.h"
#include "rna_alphabet.h"
#include "rnaforester_options.h"

// TODO in und del sind immer gleich, also kann man das in private indel fknen zusammenfassen
// TODO mit einer ordentlichen hierarchie koennte man z.b. den score in der oberklasse setzen
// man koennte dann auch einfach alle klassen in diesem file beerben, und um die affinen funktionen erweitern
// TODO zerlegen in verschiedene files, fuer profile und szalgebra extra

// Definitions and typedefs

const double DBL_NEG = -100000000.0;		// TODO the values from limits.h caused problems..
const double DBL_POS = 100000000.0;

// Class Score

// Class Score reads scoring parameters from RNAforester command line.
class Score {
	friend std::ostream& operator<< (std::ostream &out, const Score &score);
private:
    bool isAffine_;
    bool isLocal_;
    bool isDistance_;
    bool isRIBOSUM_;

public:
    double bp_rep_score_;
    double bp_indel_score_;
    double bp_indel_open_score_;
    double b_match_score_;
    double b_rep_score_;
    double b_indel_score_;
    double b_indel_open_score_;

    Score(Options &options);
    Score(const Score &s);
};


// RNA Algebra Classes


// Similarity algebra for RNA forests 

class DoubleSimiRNA_Algebra : public RNA_Algebra<double,RNA_Alphabet> {
private:
    Score s_;

public:
    double empty() const {
        return 0;
    };

    double replacepair(RNA_Alphabet la, RNA_Alphabet lb, double down, RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
        return s_.bp_rep_score_+down+over;
    };

    double replace(RNA_Alphabet a,double down, RNA_Alphabet b, double over) const {
        if (a==ALPHA_BASEPAIR && b == ALPHA_BASEPAIR)
            return s_.bp_rep_score_+down+over;
        else {
            if (a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
                return INT_MIN;
            else {
                if (a==b)
                    return s_.b_match_score_+down+over;
                else
                    return s_.b_rep_score_+down+over;
            }
        }
    };

    double del(RNA_Alphabet a,double down, double over) const {
        if (a==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double insert(double down,RNA_Alphabet b,double over) const {
        if (b==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double choice(double a, double  b) const {
        return std::max(a,b);
    };

    double worst_score() const {
        return INT_MIN;
    };

    DoubleSimiRNA_Algebra(const Score &s)
            : s_(s) {};

		// Pair indels

		// three/four possibilities for scores, choose via choice function:
		// delete pairing only, but leave bases
		// del p + rep l + mdown + rep r + over
		// delete pairing and both bases
		// del p + del l + mdown + del r + over
		// delete pairing and one base
		// del p + del l + mdown + rep r + over
		// del p + rep l + mdown + del r + over
		double deletePairOnly(RNA_Alphabet la, RNA_Alphabet lb, double mdown,   RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_rep_score_ + over;
				double result = s_.bp_indel_score_;
				if (la == lb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += mdown;
				if (ra == rb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += over;
				return result;
		}
		double deletePairAndBases(RNA_Alphabet la, double mdown, RNA_Alphabet   ra, double over) const {
				// return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}
		double deletePairAndLeftBase(RNA_Alphabet la, double mdown,                   RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_rep_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				if (ra == rb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += over;
				return result;
		}
		double deletePairAndRightBase(RNA_Alphabet la, RNA_Alphabet lb, double  mdown, RNA_Alphabet ra, double over) const {
				// return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_indel_score_ + over;
			double result = s_.bp_indel_score_;
			if (la == lb)
				result += s_.b_match_score_;
			else
				result += s_.b_rep_score_;
			result += mdown;
			result += s_.b_indel_score_;
			result += over;
			return result;
		}


		// three/four possibilities for scores, choose via choice function:
		// insert pairing only, but leave bases
		// ins p + rep l + mdown + rep r + over
		// insert pairing and both bases
		// ins p + ins l + mdown + ins r + over
		// insert pairing and one base
		// ins p + ins l + mdown + rep r + over
		// ins p + rep l + mdown + ins r + over
		double insertPairOnly(RNA_Alphabet la, RNA_Alphabet lb, double mdown,   RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_rep_score_ + over;  
				double result = s_.b_indel_score_;
				if (la == lb) {
					result += s_.b_match_score_;
				}
				else {
					result += s_.b_rep_score_;
				}
				result += mdown;
				if (ra == rb) {
					result += s_.b_match_score_;
				}
				else {
					result += s_.b_rep_score_;
				}
				result += over;
				return result;
		}
		double insertPairAndBases(RNA_Alphabet lb, double mdown, RNA_Alphabet   rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}
		double insertPairAndLeftBase(RNA_Alphabet lb, double mdown,                   RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_rep_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				if (ra == rb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += over;
				return result;
		}
		double insertPairAndRightBase(RNA_Alphabet la, RNA_Alphabet lb, double  mdown, RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				if (la == lb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}

};

class AffineDoubleSimiRNA_Algebra : public RNA_AlgebraAffine<double,RNA_Alphabet> {

private:
    Score s_;
public:
    double empty() const {
        return 0;
    };
    double replacepair(RNA_Alphabet la, RNA_Alphabet lb, double down, RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
        return s_.bp_rep_score_+down+over;
    };

    double replace(RNA_Alphabet a,double down, RNA_Alphabet b, double over) const {
        if (a==ALPHA_BASEPAIR && b == ALPHA_BASEPAIR)
            return s_.bp_rep_score_+down+over;
        else {
            if (a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
                return INT_MIN;
            else {
                if (a==b)
                    return s_.b_match_score_+down+over;
                else
                    return s_.b_rep_score_+down+over;
            }
        }
    };

    double del(RNA_Alphabet a,double down, double over) const {
        if (a==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double insert(double down,RNA_Alphabet b,double over) const {
        if (b==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double choice(double a, double  b) const {
        return std::max(a,b);
    };

    double worst_score() const {
        return INT_MIN;
    };


    double delO(RNA_Alphabet a,double down, double over) const {
        if (a==ALPHA_BASEPAIR)
            return s_.bp_indel_open_score_+down+over;
        else
            return s_.b_indel_open_score_+down+over;
    };

    double insertO(double down,RNA_Alphabet b,double over) const {
        if (b==ALPHA_BASEPAIR)
            return s_.bp_indel_open_score_+down+over;
        else
            return s_.b_indel_open_score_+down+over;
    };

    AffineDoubleSimiRNA_Algebra(const Score &s)
            : s_(s) {};

		// Pair indels

									
		// three/four possibilities for scores, choose via choice function:
		// delete pairing only, but leave bases
		// del p + rep l + mdown + rep r + over
		// delete pairing and both bases
		// del p + del l + mdown + del r + over
		// delete pairing and one base
		// del p + del l + mdown + rep r + over
		// del p + rep l + mdown + del r + over
		double deletePairOnly(RNA_Alphabet la, RNA_Alphabet lb, double    mdown,
												RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
				double result = s_.bp_indel_score_;
				if (la == lb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += mdown;
				if (ra == rb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += over;
				return result;
				// return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_rep_score_ + over;
		}

		double deletePairAndBases(RNA_Alphabet la, double mdown,
												RNA_Alphabet ra, double over) const {
				// return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}

		double deletePairAndLeftBase(RNA_Alphabet la, double mdown,
												RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_rep_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				if (ra == rb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += over;
				return result;

		}

		double deletePairAndRightBase(RNA_Alphabet la, RNA_Alphabet lb,   double mdown,
												RNA_Alphabet ra, double over) const {
				// return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_indel_score_ + over;
			double result = s_.bp_indel_score_;
			if (la == lb)
				result += s_.b_match_score_;
			else
				result += s_.b_rep_score_;
			result += mdown;
			result += s_.b_indel_score_;
			result += over;
			return result;
		}


		// three/four possibilities for scores, choose via choice function:
		// insert pairing only, but leave bases
		// ins p + rep l + mdown + rep r + over
		// insert pairing and both bases
		// ins p + ins l + mdown + ins r + over
		// insert pairing and one base
		// ins p + ins l + mdown + rep r + over
		// ins p + rep l + mdown + ins r + over
		double insertPairOnly(RNA_Alphabet la, RNA_Alphabet lb, double    mdown,
												RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_rep_score_ + over; 
				double result = s_.b_indel_score_;
				if (la == lb) {
					result += s_.b_match_score_;
				}
				else {
					result += s_.b_rep_score_;
				}
				result += mdown;
				if (ra == rb) {
					result += s_.b_match_score_;
				}
				else {
					result += s_.b_rep_score_;
				}
				result += over;
				return result;
		}
		double insertPairAndBases(RNA_Alphabet lb, double mdown,
												RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}
		double insertPairAndLeftBase(RNA_Alphabet lb, double mdown,
												RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_rep_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				if (ra == rb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += over;
				return result;
		}
		double insertPairAndRightBase(RNA_Alphabet la, RNA_Alphabet lb,   double mdown,
												RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				if (la == lb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}

};

// Distance algebra for RNA forests 
class DoubleDistRNA_Algebra : public RNA_Algebra<double,RNA_Alphabet> {
private:
    Score s_;

public:
    double empty() const {
        return 0;
    };
    double replacepair(RNA_Alphabet la, RNA_Alphabet lb, double down, RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
        return s_.bp_rep_score_+down+over;
    };

    double replace(RNA_Alphabet a,double down, RNA_Alphabet b, double over) const {
        if (a==ALPHA_BASEPAIR && b == ALPHA_BASEPAIR)
            return s_.bp_rep_score_+down+over;
        else {
            if (a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
                return INT_MAX;
            else {
                if (a==b)
                    return s_.b_match_score_+down+over;
                else
                    return s_.b_rep_score_+down+over;
            }
        }
    };

    double del(RNA_Alphabet a,double down,double over) const {
        if (a==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double insert(double down,RNA_Alphabet b,double over) const {
        if (b==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double choice(double a, double  b) const {
        return std::min(a,b);
    };

    double worst_score() const {
        return INT_MAX;
    };

    DoubleDistRNA_Algebra(const Score &s)
            : s_(s) {};
						
		// Pair indel
		// three/four possibilities for scores, choose via choice function:
		// delete pairing only, but leave bases
		// del p + rep l + mdown + rep r + over
		// delete pairing and both bases
		// del p + del l + mdown + del r + over
		// delete pairing and one base
		// del p + del l + mdown + rep r + over
		// del p + rep l + mdown + del r + over
		double deletePairOnly(RNA_Alphabet la, RNA_Alphabet lb, double mdown,   RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
// 				return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_rep_score_ + over;
				double result = s_.bp_indel_score_;
				if (la == lb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += mdown;
				if (ra == rb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += over;
				return result;
		}
		double deletePairAndBases(RNA_Alphabet la, double mdown, RNA_Alphabet   ra, double over) const {
// 				return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}
		double deletePairAndLeftBase(RNA_Alphabet la, double mdown,                   RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
// 				return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_rep_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				if (ra == rb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += over;
				return result;
		}
		double deletePairAndRightBase(RNA_Alphabet la, RNA_Alphabet lb, double  mdown, RNA_Alphabet ra, double over) const {
// 				return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_indel_score_ + over;
			double result = s_.bp_indel_score_;
			if (la == lb)
				result += s_.b_match_score_;
			else
				result += s_.b_rep_score_;
			result += mdown;
			result += s_.b_indel_score_;
			result += over;
			return result;

		}

		// three/four possibilities for scores, choose via choice function:
		// insert pairing only, but leave bases
		// ins p + rep l + mdown + rep r + over
		// insert pairing and both bases
		// ins p + ins l + mdown + ins r + over
		// insert pairing and one base
		// ins p + ins l + mdown + rep r + over
		// ins p + rep l + mdown + ins r + over
		double insertPairOnly(RNA_Alphabet la, RNA_Alphabet lb, double mdown,   RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
				return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_rep_score_ + over; 
				double result = s_.b_indel_score_;
				if (la == lb) {
					result += s_.b_match_score_;
				}
				else {
					result += s_.b_rep_score_;
				}
				result += mdown;
				if (ra == rb) {
					result += s_.b_match_score_;
				}
				else {
					result += s_.b_rep_score_;
				}
				result += over;
				return result;
		}
		double insertPairAndBases(RNA_Alphabet lb, double mdown, RNA_Alphabet   rb, double over) const {
// 				return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}
		double insertPairAndLeftBase(RNA_Alphabet lb, double mdown,                   RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
// 				return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_rep_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				if (ra == rb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += over;
				return result;
		}
		double insertPairAndRightBase(RNA_Alphabet la, RNA_Alphabet lb, double  mdown, RNA_Alphabet rb, double over) const {
// 				return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				if (la == lb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}


};

// Distance algebra for RNA forests 
class AffineDoubleDistRNA_Algebra : public RNA_AlgebraAffine<double,RNA_Alphabet>  {
private:
    Score s_;
public:
    double empty() const {
        return 0;
    };
    double replacepair(RNA_Alphabet la, RNA_Alphabet lb, double down, RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
        return s_.bp_rep_score_+down+over;
    };

    double replace(RNA_Alphabet a,double down, RNA_Alphabet b, double over) const {
        if (a==ALPHA_BASEPAIR && b == ALPHA_BASEPAIR)
            return s_.bp_rep_score_+down+over;
        else {
            if (a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
                return INT_MAX;
            else {
                if (a==b)
                    return s_.b_match_score_+down+over;
                else
                    return s_.b_rep_score_+down+over;
            }
        }
    };

    double del(RNA_Alphabet a,double down,double over) const {
        if (a==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double insert(double down,RNA_Alphabet b,double over) const {
        if (b==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double choice(double a, double  b) const {
        return std::min(a,b);
    };

    double worst_score() const {
        return INT_MAX;
    };


    double delO(RNA_Alphabet a,double down,double over) const {
        if (a==ALPHA_BASEPAIR)
            return s_.bp_indel_open_score_+down+over;
        else
            return s_.b_indel_open_score_+down+over;
    };

    double insertO(double down,RNA_Alphabet b,double over) const {
        if (b==ALPHA_BASEPAIR)
            return s_.bp_indel_open_score_+down+over;
        else
            return s_.b_indel_open_score_+down+over;
    };

    AffineDoubleDistRNA_Algebra(const Score &s)
            : s_(s) {};

		// three/four possibilities for scores, choose via choice function:
		// delete pairing only, but leave bases
		// del p + rep l + mdown + rep r + over
		// delete pairing and both bases
		// del p + del l + mdown + del r + over
		// delete pairing and one base
		// del p + del l + mdown + rep r + over
		// del p + rep l + mdown + del r + over
		double deletePairOnly(RNA_Alphabet la, RNA_Alphabet lb, double    mdown,
										RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_rep_score_ + over;
				double result = s_.bp_indel_score_;
				if (la == lb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += mdown;
				if (ra == rb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += over;
				return result;
		}
		double deletePairAndBases(RNA_Alphabet la, double mdown,
										RNA_Alphabet ra, double over) const {
				// return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}
		double deletePairAndLeftBase(RNA_Alphabet la, double mdown,
										RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_rep_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				if (ra == rb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += over;
				return result;
		}
		double deletePairAndRightBase(RNA_Alphabet la, RNA_Alphabet lb,   double mdown,
										RNA_Alphabet ra, double over) const {
				// return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_indel_score_ + over;
			double result = s_.bp_indel_score_;
			if (la == lb)
				result += s_.b_match_score_;
			else
				result += s_.b_rep_score_;
			result += mdown;
			result += s_.b_indel_score_;
			result += over;
			return result;
		}

		// three/four possibilities for scores, choose via choice function:
		// insert pairing only, but leave bases
		// ins p + rep l + mdown + rep r + over
		// insert pairing and both bases
		// ins p + ins l + mdown + ins r + over
		// insert pairing and one base
		// ins p + ins l + mdown + rep r + over
		// ins p + rep l + mdown + ins r + over
		double insertPairOnly(RNA_Alphabet la, RNA_Alphabet lb, double    mdown,
												RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_rep_score_ + over; 
				double result = s_.b_indel_score_;
				if (la == lb) {
					result += s_.b_match_score_;
				}
				else {
					result += s_.b_rep_score_;
				}
				result += mdown;
				if (ra == rb) {
					result += s_.b_match_score_;
				}
				else {
					result += s_.b_rep_score_;
				}
				result += over;
				return result;
		}
		double insertPairAndBases(RNA_Alphabet lb, double mdown,
												RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}
		double insertPairAndLeftBase(RNA_Alphabet lb, double mdown,
												RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_rep_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				if (ra == rb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += over;
				return result;
		}
		double insertPairAndRightBase(RNA_Alphabet la, RNA_Alphabet lb,   double mdown,
												RNA_Alphabet rb, double over) const {
				// return s_.bp_indel_score_ + s_.b_rep_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				if (la == lb)
					result += s_.b_match_score_;
				else
					result += s_.b_rep_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}


};

/** RIBOSUM85-60 matrix published in RSEARCH: Finding homologs of single structured RNA sequences
 *  R. Klein and S. Eddy, BMC Bioinformatics 2003 Vol.4
 */
class RIBOSUM8560 : public RNA_Algebra<double,RNA_Alphabet> {
private:
    double baseSubstMtrx_[4][4];
    double basepairSubstMtrx_[4][4][4][4];
    Score s_;
public:
    double empty() const {
        return 0;
    };
    double replacepair(RNA_Alphabet la, RNA_Alphabet lb, double down, RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
        int i,j,k,l;
        i = alpha2RNA_Alpha(la);
        j = alpha2RNA_Alpha(ra);
        k = alpha2RNA_Alpha(lb);
        l = alpha2RNA_Alpha(rb);

        return basepairSubstMtrx_[i][j][k][l]+down+over;
    };

    double replace(RNA_Alphabet a,double down, RNA_Alphabet b, double over) const {
        assert(!(a==ALPHA_BASEPAIR && b==ALPHA_BASEPAIR));

        if (a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
            return INT_MIN;
        else {
            int i,j;
            i = alpha2RNA_Alpha(a);
            j = alpha2RNA_Alpha(b);

            return baseSubstMtrx_[i][j]+down+over;
        }
    };

    double del(RNA_Alphabet a,double down, double over) const {
        if (a==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double insert(double down,RNA_Alphabet b,double over) const {
        if (b==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double choice(double a, double  b) const {
        return std::max(a,b);
    };

    double worst_score() const {
        return INT_MIN;
    };

    RIBOSUM8560(const Score &s)
            : s_(s) {
        int i,j,k,l;

        // set substitution matrices

        // base replacement
        baseSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_A]=222;
        baseSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_C]=-186;
        baseSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_G]=-146;
        baseSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U]=-139;

        baseSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_C]=116;
        baseSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G]=-248;
        baseSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_U]=-105;

        baseSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_G]=103;
        baseSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=-174;

        baseSubstMtrx_[ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U]=165;

        // copy triangle
        for (i=0; i<=ALPHA_PRO_BASE_U; i++)
            for (j=0; j<i; j++)
                baseSubstMtrx_[i][j]=baseSubstMtrx_[j][i];

        // basepair replacement

        // set default score. This score should never be used since the scores for canonical basepairs are defined later
        for (i=0; i<=ALPHA_PRO_BASE_U; i++)
            for (j=0; j<=ALPHA_PRO_BASE_U; j++)
                for (k=i; k<=ALPHA_PRO_BASE_U; k++)
                    for (l=j; l<=ALPHA_PRO_BASE_U; l++)
                        basepairSubstMtrx_[i][j][k][l]=-1000;

        basepairSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U]=449;
        basepairSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G]=167;
        basepairSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C]=270;
        basepairSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=59;
        basepairSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=161;
        basepairSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=-51;

        basepairSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G]=536;
        basepairSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C]=211;
        basepairSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=-27;
        basepairSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=275;
        basepairSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=132;

        basepairSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C]=562;
        basepairSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=121;
        basepairSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=167;
        basepairSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=-8;

        basepairSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=347;
        basepairSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=-57;
        basepairSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=-209;

        basepairSubstMtrx_[ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=497;
        basepairSubstMtrx_[ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=114;

        basepairSubstMtrx_[ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=336;

        // copy triangle
        for (i=0; i<=ALPHA_PRO_BASE_U; i++)
            for (j=0; j<=ALPHA_PRO_BASE_U; j++)
                for (k=0; k<=ALPHA_PRO_BASE_U; k++)
                    for (l=0; l<=ALPHA_PRO_BASE_U; l++)
                        if (k<i || (k==i && l<j))
                            basepairSubstMtrx_[i][j][k][l] = basepairSubstMtrx_[k][l][i][j];
    };

		// Pair indel
						
		// three/four possibilities for scores, choose via choice function:
		// delete pairing only, but leave bases
		// del p + rep l + mdown + rep r + over
		// delete pairing and both bases
		// del p + del l + mdown + del r + over
		// delete pairing and one base
		// del p + del l + mdown + rep r + over
		// del p + rep l + mdown + del r + over
		double deletePairOnly(const RNA_Alphabet la, const RNA_Alphabet lb, const double mdown,
										const RNA_Alphabet ra, const RNA_Alphabet rb, const double over) const {

				int i, j, k, l;
				i = alpha2RNA_Alpha(la);
				j = alpha2RNA_Alpha(ra);
				k = alpha2RNA_Alpha(lb);
				l = alpha2RNA_Alpha(rb);

				return s_.bp_indel_score_ + basepairSubstMtrx_[i][j][k][l] + mdown + over;
		}
		double deletePairAndBases(const RNA_Alphabet la, const double mdown,
										const RNA_Alphabet ra, const double over) const {

// 				return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}
		double deletePairAndLeftBase(const RNA_Alphabet la, const double mdown,
										const RNA_Alphabet ra, const RNA_Alphabet rb, const double over) const {

				int  k, l;
				k = alpha2RNA_Alpha(ra);
				l = alpha2RNA_Alpha(rb);

				return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + baseSubstMtrx_[k][l] + over;
		}
		double deletePairAndRightBase(const RNA_Alphabet la, const RNA_Alphabet lb, const double mdown,
										const RNA_Alphabet ra, const double over) const {

				int i, j;
				i = alpha2RNA_Alpha(la);
				j = alpha2RNA_Alpha(ra);

				return s_.bp_indel_score_ + baseSubstMtrx_[i][j] + mdown + s_.b_indel_score_ + over;
		}

		// three/four possibilities for scores, choose via choice function:
		// insert pairing only, but leave bases
		// ins p + rep l + mdown + rep r + over
		// insert pairing and both bases
		// ins p + ins l + mdown + ins r + over
		// insert pairing and one base
		// ins p + ins l + mdown + rep r + over
		// ins p + rep l + mdown + ins r + over
		double insertPairOnly(const RNA_Alphabet la, const RNA_Alphabet lb, const double mdown,
										const RNA_Alphabet ra, const RNA_Alphabet rb, const double over) const {

				int i, j, k, l;
				i = alpha2RNA_Alpha(la);
				j = alpha2RNA_Alpha(ra);
				k = alpha2RNA_Alpha(lb);
				l = alpha2RNA_Alpha(rb);

				return s_.bp_indel_score_ + basepairSubstMtrx_[i][j][k][l] + mdown + over;
		}

		double insertPairAndBases(const RNA_Alphabet lb, const double mdown,
										const RNA_Alphabet rb, const double over) const {

// 				return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}

		double insertPairAndLeftBase(const RNA_Alphabet lb, const double mdown,
										const RNA_Alphabet ra, const RNA_Alphabet rb, const double over) const {

				int k, l;
				k = alpha2RNA_Alpha(lb);
				l = alpha2RNA_Alpha(rb);

				return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + baseSubstMtrx_[k][l] + over;
		}

		double insertPairAndRightBase(const RNA_Alphabet la, const RNA_Alphabet lb, const double mdown,
										const RNA_Alphabet rb, const double over) const {

				int i, j;
				i = alpha2RNA_Alpha(lb);
				j = alpha2RNA_Alpha(rb);

				return s_.bp_indel_score_ + baseSubstMtrx_[i][j] + mdown + s_.b_indel_score_ + over;
		}

	
};

class AffineRIBOSUM8560 : public RNA_AlgebraAffine<double,RNA_Alphabet> {
private:
    Score s_;
    double baseSubstMtrx_[4][4];
    double basepairSubstMtrx_[4][4][4][4];
public:
    double empty() const {
        return 0;
    };
    double replacepair(RNA_Alphabet la, RNA_Alphabet lb, double down, RNA_Alphabet ra, RNA_Alphabet rb, double over) const {
        int i,j,k,l;
        i = alpha2RNA_Alpha(la);
        j = alpha2RNA_Alpha(ra);
        k = alpha2RNA_Alpha(lb);
        l = alpha2RNA_Alpha(rb);

        return basepairSubstMtrx_[i][j][k][l]+down+over;
    };

    double replace(RNA_Alphabet a,double down, RNA_Alphabet b, double over) const {
        assert(!(a==ALPHA_BASEPAIR && b==ALPHA_BASEPAIR));

        if (a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
            return INT_MIN;
        else {
            int i,j;
            i = alpha2RNA_Alpha(a);
            j = alpha2RNA_Alpha(b);

            return baseSubstMtrx_[i][j]+down+over;
        }
    };

    double del(RNA_Alphabet a,double down, double over) const {
        if (a==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double insert(double down,RNA_Alphabet b,double over) const {
        if (b==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double choice(double a, double  b) const {
        return std::max(a,b);
    };

    double worst_score() const {
        return INT_MIN;
    };

    AffineRIBOSUM8560(const Score &s)
            : s_(s) {
        int i,j,k,l;

        // set substitution matrices

        // base replacement
        baseSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_A]=222;
        baseSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_C]=-186;
        baseSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_G]=-146;
        baseSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U]=-139;

        baseSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_C]=116;
        baseSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G]=-248;
        baseSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_U]=-105;

        baseSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_G]=103;
        baseSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=-174;

        baseSubstMtrx_[ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U]=165;

        // copy triangle
        for (i=0; i<=ALPHA_PRO_BASE_U; i++)
            for (j=0; j<i; j++)
                baseSubstMtrx_[i][j]=baseSubstMtrx_[j][i];

        // basepair replacement

        // set default score. This score should never be used since the scores for canonical basepairs are defined later
        for (i=0; i<=ALPHA_PRO_BASE_U; i++)
            for (j=0; j<=ALPHA_PRO_BASE_U; j++)
                for (k=i; k<=ALPHA_PRO_BASE_U; k++)
                    for (l=j; l<=ALPHA_PRO_BASE_U; l++)
                        basepairSubstMtrx_[i][j][k][l]=-1000;

        basepairSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U]=449;
        basepairSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G]=167;
        basepairSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C]=270;
        basepairSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=59;
        basepairSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=161;
        basepairSubstMtrx_[ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=-51;

        basepairSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G]=536;
        basepairSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C]=211;
        basepairSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=-27;
        basepairSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=275;
        basepairSubstMtrx_[ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=132;

        basepairSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C]=562;
        basepairSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=121;
        basepairSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=167;
        basepairSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_C][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=-8;

        basepairSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U]=347;
        basepairSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=-57;
        basepairSubstMtrx_[ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=-209;

        basepairSubstMtrx_[ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A]=497;
        basepairSubstMtrx_[ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_A][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=114;

        basepairSubstMtrx_[ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G][ALPHA_PRO_BASE_U][ALPHA_PRO_BASE_G]=336;

        // copy triangle
        for (i=0; i<=ALPHA_PRO_BASE_U; i++)
            for (j=0; j<=ALPHA_PRO_BASE_U; j++)
                for (k=0; k<=ALPHA_PRO_BASE_U; k++)
                    for (l=0; l<=ALPHA_PRO_BASE_U; l++)
                        if (k<i || (k==i && l<j))
                            basepairSubstMtrx_[i][j][k][l]=basepairSubstMtrx_[k][l][i][j];
    };

    double delO(RNA_Alphabet a,double down, double over) const {
        if (a==ALPHA_BASEPAIR)
            return s_.bp_indel_open_score_+down+over;
        else
            return s_.b_indel_open_score_+down+over;
    };

    double insertO(double down,RNA_Alphabet b,double over) const {
        if (b==ALPHA_BASEPAIR)
            return s_.bp_indel_open_score_+down+over;
        else
            return s_.b_indel_open_score_+down+over;
    };

		// Pair indel

		// three/four possibilities for scores, choose via choice function:
		// delete pairing only, but leave bases
		// del p + rep l + mdown + rep r + over
		// delete pairing and both bases
		// del p + del l + mdown + del r + over
		// delete pairing and one base
		// del p + del l + mdown + rep r + over
		// del p + rep l + mdown + del r + over
		double deletePairOnly(const RNA_Alphabet la, const RNA_Alphabet lb, const double mdown,
										const RNA_Alphabet ra, const RNA_Alphabet rb, const double over) const {

				int i, j, k, l;
				i = alpha2RNA_Alpha(la);
				j = alpha2RNA_Alpha(ra);
				k = alpha2RNA_Alpha(lb);
				l = alpha2RNA_Alpha(rb);

				return s_.bp_indel_score_ + basepairSubstMtrx_[i][j][k][l] + mdown + over;
		}
		double deletePairAndBases(const RNA_Alphabet la, const double mdown,
										const RNA_Alphabet ra, const double over) const {

// 				return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}
		double deletePairAndLeftBase(const RNA_Alphabet la, const double mdown,
										const RNA_Alphabet ra, const RNA_Alphabet rb, const double over) const {

				int  k, l;
				k = alpha2RNA_Alpha(ra);
				l = alpha2RNA_Alpha(rb);

				return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + baseSubstMtrx_[k][l] + over;
		}
		double deletePairAndRightBase(const RNA_Alphabet la, const RNA_Alphabet lb, const double mdown,
										const RNA_Alphabet ra, const double over) const {

				int i, j;
				i = alpha2RNA_Alpha(la);
				j = alpha2RNA_Alpha(ra);

				return s_.bp_indel_score_ + baseSubstMtrx_[i][j] + mdown + s_.b_indel_score_ + over;
		}

		// three/four possibilities for scores, choose via choice function:
		// insert pairing only, but leave bases
		// ins p + rep l + mdown + rep r + over
		// insert pairing and both bases
		// ins p + ins l + mdown + ins r + over
		// insert pairing and one base
		// ins p + ins l + mdown + rep r + over
		// ins p + rep l + mdown + ins r + over
		double insertPairOnly(const RNA_Alphabet la, const RNA_Alphabet lb, const double mdown,
										const RNA_Alphabet ra, const RNA_Alphabet rb, const double over) const {

				int i, j, k, l;
				i = alpha2RNA_Alpha(la);
				j = alpha2RNA_Alpha(ra);
				k = alpha2RNA_Alpha(lb);
				l = alpha2RNA_Alpha(rb);

				return s_.bp_indel_score_ + basepairSubstMtrx_[i][j][k][l] + mdown + over;
		}

		double insertPairAndBases(const RNA_Alphabet lb, const double mdown,
										const RNA_Alphabet rb, const double over) const {

// 				return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + s_.b_indel_score_ + over;
				double result = s_.bp_indel_score_;
				result += s_.b_indel_score_;
				result += mdown;
				result += s_.b_indel_score_;
				result += over;
				return result;
		}

		double insertPairAndLeftBase(const RNA_Alphabet lb, const double mdown,
										const RNA_Alphabet ra, const RNA_Alphabet rb, const double over) const {

				int k, l;
				k = alpha2RNA_Alpha(lb);
				l = alpha2RNA_Alpha(rb);

				return s_.bp_indel_score_ + s_.b_indel_score_ + mdown + baseSubstMtrx_[k][l] + over;
		}

		double insertPairAndRightBase(const RNA_Alphabet la, const RNA_Alphabet lb, const double mdown,
										const RNA_Alphabet rb, const double over) const {

				int i, j;
				i = alpha2RNA_Alpha(lb);
				j = alpha2RNA_Alpha(rb);

				return s_.bp_indel_score_ + baseSubstMtrx_[i][j] + mdown + s_.b_indel_score_ + over;
		}

};




// RNA Profile Algebra Classes

// Similarity algebra for RNA profile forests
class DoubleSimiProfileAlgebra : public Algebra<double,RNA_Alphabet_Profile> {
private:
    Score s_;

public:
    double empty() const {
        return 0.0;
    };
    double replace(RNA_Alphabet_Profile a,double down, RNA_Alphabet_Profile b, double over) const {
        if (a.p[ALPHA_PRO_BASEPAIR]>0 && b.p[ALPHA_PRO_BASEPAIR]>0) {
            // pair replacement
            return a.p[ALPHA_PRO_BASEPAIR]*b.p[ALPHA_PRO_BASEPAIR]*s_.bp_rep_score_ +
                   down+over;
        } else {
            if (a.p[ALPHA_PRO_BASE]>0 && b.p[ALPHA_PRO_BASE]>0) {
                double s=0;

                // base replacement
                for (int i=ALPHA_PRO_BASE_A; i<=ALPHA_PRO_BASE_U; i++)
                    for (int j=ALPHA_PRO_BASE_A; j<=ALPHA_PRO_BASE_U; j++)
                        s+= i==j ? a.p[i]*b.p[j]*s_.b_match_score_ : a.p[i]*b.p[j]*s_.b_rep_score_;

                if (s==0) // no sequence information
                    s=a.p[ALPHA_PRO_BASE]*b.p[ALPHA_PRO_BASE]*s_.b_rep_score_;

                return s+down+over;
            } else {
                // undefined operation (replace base by basepair ??)
                return DBL_NEG/4;
            }
        }
    };

    double del(RNA_Alphabet_Profile a,double down, double over) const {
        if (a.p[ALPHA_PRO_BASEPAIR]>0)
            return a.p[ALPHA_PRO_BASEPAIR]*s_.bp_indel_score_+down+over;
        else
            return a.p[ALPHA_PRO_BASE]*s_.b_indel_score_+down+over;
    };

    double insert(double down,RNA_Alphabet_Profile b,double over) const {
        if (b.p[ALPHA_PRO_BASEPAIR]>0)
            return b.p[ALPHA_PRO_BASEPAIR]*s_.bp_indel_score_+down+over;
        else
            return b.p[ALPHA_PRO_BASE]*s_.b_indel_score_+down+over;
    };

    double choice(double a, double  b) const {
        return std::max(a,b);
    };

    double worst_score() const {
        return DBL_NEG;
    };

    DoubleSimiProfileAlgebra(const Score &s)
            : s_(s) {};
};

class AffineDoubleSimiProfileAlgebra : public AlgebraAffine<double,RNA_Alphabet_Profile> {
private:
    Score s_;
public:
    double empty() const {
        return 0.0;
    };
    double replace(RNA_Alphabet_Profile a,double down, RNA_Alphabet_Profile b, double over) const {
        if (a.p[ALPHA_PRO_BASEPAIR]>0 && b.p[ALPHA_PRO_BASEPAIR]>0) {
            // pair replacement
            return a.p[ALPHA_PRO_BASEPAIR]*b.p[ALPHA_PRO_BASEPAIR]*s_.bp_rep_score_ +
                   down+over;
        } else {
            if (a.p[ALPHA_PRO_BASE]>0 && b.p[ALPHA_PRO_BASE]>0) {
                double s=0;

                // base replacement
                for (int i=ALPHA_PRO_BASE_A; i<=ALPHA_PRO_BASE_U; i++)
                    for (int j=ALPHA_PRO_BASE_A; j<=ALPHA_PRO_BASE_U; j++)
                        s+= i==j ? a.p[i]*b.p[j]*s_.b_match_score_ : a.p[i]*b.p[j]*s_.b_rep_score_;

                if (s==0) // no sequence information
                    s=a.p[ALPHA_PRO_BASE]*b.p[ALPHA_PRO_BASE]*s_.b_rep_score_;

                return s+down+over;
            } else {
                // undefined operation (replace base by basepair ??)
                return DBL_NEG/4;
            }
        }
    };

    double del(RNA_Alphabet_Profile a,double down, double over) const {
        if (a.p[ALPHA_PRO_BASEPAIR]>0)
            return a.p[ALPHA_PRO_BASEPAIR]*s_.bp_indel_score_+down+over;
        else
            return a.p[ALPHA_PRO_BASE]*s_.b_indel_score_+down+over;
    };

    double insert(double down,RNA_Alphabet_Profile b,double over) const {
        if (b.p[ALPHA_PRO_BASEPAIR]>0)
            return b.p[ALPHA_PRO_BASEPAIR]*s_.bp_indel_score_+down+over;
        else
            return b.p[ALPHA_PRO_BASE]*s_.b_indel_score_+down+over;
    };

    double choice(double a, double  b) const {
        return std::max(a,b);
    };

    double worst_score() const {
        return DBL_NEG;
    };


    double delO(RNA_Alphabet_Profile a,double down, double over) const {
        if (a.p[ALPHA_PRO_BASEPAIR]>0)
            return a.p[ALPHA_PRO_BASEPAIR]*s_.bp_indel_open_score_+down+over;
        else
            return a.p[ALPHA_PRO_BASE]*s_.b_indel_open_score_+down+over;
    };

    double insertO(double down,RNA_Alphabet_Profile b,double over) const {
        if (b.p[ALPHA_PRO_BASEPAIR]>0)
            return b.p[ALPHA_PRO_BASEPAIR]*s_.bp_indel_open_score_+down+over;
        else
            return b.p[ALPHA_PRO_BASE]*s_.b_indel_open_score_+down+over;
    };

    AffineDoubleSimiProfileAlgebra(const Score &s)
            : s_(s) {};
};




// Distance algebra for RNA profile forests
class DoubleDistProfileAlgebra : public Algebra<double,RNA_Alphabet_Profile> {
private:
    Score s_;

public:
    double empty() const {
        return 0.0;
    };
    double replace(RNA_Alphabet_Profile a,double down, RNA_Alphabet_Profile b, double over) const {
        TRACE(DBG_ALGEBRA,"rep","inside!!!");

        if (a.p[ALPHA_PRO_BASEPAIR]>0 && b.p[ALPHA_PRO_BASEPAIR]>0) {
            // pair replacement
            return a.p[ALPHA_PRO_BASEPAIR]*b.p[ALPHA_PRO_BASEPAIR]*s_.bp_rep_score_ +
                   down+over;
        } else {
            if (a.p[ALPHA_PRO_BASE]>0 && b.p[ALPHA_PRO_BASE]>0) {
                double s=0;

                // base replacement
                for (int i=ALPHA_PRO_BASE_A; i<=ALPHA_PRO_BASE_U; i++)
                    for (int j=ALPHA_PRO_BASE_A; j<=ALPHA_PRO_BASE_U; j++)
                        s+= i==j ? a.p[i]*b.p[j]*s_.b_match_score_ : a.p[i]*b.p[j]*s_.b_rep_score_;

                if (s==0) // no sequence information
                    s=a.p[ALPHA_PRO_BASE]*b.p[ALPHA_PRO_BASE]*s_.b_rep_score_;

                return s+down+over;
            } else {
                // undefined operation (replace base by basepair ??)
                return DBL_POS/4;
            }
        }
    };

    double del(RNA_Alphabet_Profile a,double down, double over) const {
        if (a.p[ALPHA_PRO_BASEPAIR]>0)
            return a.p[ALPHA_PRO_BASEPAIR]*s_.bp_indel_score_+down+over;
        else
            return a.p[ALPHA_PRO_BASE]*s_.b_indel_score_+down+over;
    };

    double insert(double down,RNA_Alphabet_Profile b,double over) const {
        if (b.p[ALPHA_PRO_BASEPAIR]>0)
            return b.p[ALPHA_PRO_BASEPAIR]*s_.bp_indel_score_+down+over;
        else
            return b.p[ALPHA_PRO_BASE]*s_.b_indel_score_+down+over;
    };

    double choice(double a, double  b) const {
        return std::min(a,b);
    };

    double worst_score() const {
        return DBL_POS;
    };

    DoubleDistProfileAlgebra(const Score &s)
            : s_(s) {};
};


// Distance algebra for RNA profile forests
class AffineDoubleDistProfileAlgebra : public AlgebraAffine<double,RNA_Alphabet_Profile> {
private:
    Score s_;
public:
    double empty() const {
        return 0.0;
    };
    double replace(RNA_Alphabet_Profile a,double down, RNA_Alphabet_Profile b, double over) const {
        TRACE(DBG_ALGEBRA,"rep","inside!!!");

        if (a.p[ALPHA_PRO_BASEPAIR]>0 && b.p[ALPHA_PRO_BASEPAIR]>0) {
            // pair replacement
            return a.p[ALPHA_PRO_BASEPAIR]*b.p[ALPHA_PRO_BASEPAIR]*s_.bp_rep_score_ +
                   down+over;
        } else {
            if (a.p[ALPHA_PRO_BASE]>0 && b.p[ALPHA_PRO_BASE]>0) {
                double s=0;

                // base replacement
                for (int i=ALPHA_PRO_BASE_A; i<=ALPHA_PRO_BASE_U; i++)
                    for (int j=ALPHA_PRO_BASE_A; j<=ALPHA_PRO_BASE_U; j++)
                        s+= i==j ? a.p[i]*b.p[j]*s_.b_match_score_ : a.p[i]*b.p[j]*s_.b_rep_score_;

                if (s==0) // no sequence information
                    s=a.p[ALPHA_PRO_BASE]*b.p[ALPHA_PRO_BASE]*s_.b_rep_score_;

                return s+down+over;
            } else {
                // undefined operation (replace base by basepair ??)
                return DBL_POS/4;
            }
        }
    };

    double del(RNA_Alphabet_Profile a,double down, double over) const {
        if (a.p[ALPHA_PRO_BASEPAIR]>0)
            return a.p[ALPHA_PRO_BASEPAIR]*s_.bp_indel_score_+down+over;
        else
            return a.p[ALPHA_PRO_BASE]*s_.b_indel_score_+down+over;
    };

    double insert(double down,RNA_Alphabet_Profile b,double over) const {
        if (b.p[ALPHA_PRO_BASEPAIR]>0)
            return b.p[ALPHA_PRO_BASEPAIR]*s_.bp_indel_score_+down+over;
        else
            return b.p[ALPHA_PRO_BASE]*s_.b_indel_score_+down+over;
    };

    double choice(double a, double  b) const {
        return std::min(a,b);
    };

    double worst_score() const {
        return DBL_POS;
    };


    double delO(RNA_Alphabet_Profile a,double down, double over) const {
        if (a.p[ALPHA_PRO_BASEPAIR]>0)
            return a.p[ALPHA_PRO_BASEPAIR]*s_.bp_indel_open_score_+down+over;
        else
            return a.p[ALPHA_PRO_BASE]*s_.b_indel_open_score_+down+over;
    };

    double insertO(double down,RNA_Alphabet_Profile b,double over) const {
        if (b.p[ALPHA_PRO_BASEPAIR]>0)
            return b.p[ALPHA_PRO_BASEPAIR]*s_.bp_indel_open_score_+down+over;
        else
            return b.p[ALPHA_PRO_BASE]*s_.b_indel_open_score_+down+over;
    };

    AffineDoubleDistProfileAlgebra(const Score &s)
            : s_(s) {};
};

// SZAlgebra Classes

class IntSimiSZAlgebra : public SZAlgebra<double,RNA_Alphabet> {
private:
    Score s_;

public:
    double empty() const {
        return 0;
    };
    double replace(RNA_Alphabet a,double down, RNA_Alphabet b) const {
        if (a==ALPHA_BASEPAIR && b == ALPHA_BASEPAIR)
            return s_.bp_rep_score_+down;
        else {
            if (a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
                return INT_MIN;
            else {
                if (a==b)
                    return s_.b_match_score_+down;
                else
                    return s_.b_rep_score_+down;
            }
        }
    };

    double del(RNA_Alphabet a,double down) const {
        if (a==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down;
        else
            return s_.b_indel_score_+down;
    };

    double insert(double down,RNA_Alphabet b) const {
        if (b==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down;
        else
            return s_.b_indel_score_+down;
    };

    double choice(double a, double  b) const {
        return std::max(a,b);
    };

    double worst_score() const {
        return INT_MIN;
    };

    IntSimiSZAlgebra(const Score &s)
            : s_(s) {};
};


class IntDistSZAlgebra : public SZAlgebra<double,RNA_Alphabet> {
private:
    Score s_;

public:
    double empty() const {
        return 0;
    };
    double replace(RNA_Alphabet a,double down, RNA_Alphabet b) const {
        if (a==ALPHA_BASEPAIR && b == ALPHA_BASEPAIR)
            return s_.bp_rep_score_+down;
        else {
            if (a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
                return INT_MAX;
            else {
                if (a==b)
                    return s_.b_match_score_+down;
                else
                    return s_.b_rep_score_+down;
            }
        }
    };

    double del(RNA_Alphabet a,double down) const {
        if (a==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down;
        else
            return s_.b_indel_score_+down;
    };

    double insert(double down,RNA_Alphabet b) const {
        if (b==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down;
        else
            return s_.b_indel_score_+down;
    };

    double choice(double a, double  b) const {
        return std::min(a,b);
    };

    double worst_score() const {
        return INT_MAX;
    };

    IntDistSZAlgebra(const Score &s)
            : s_(s) {};
};

// General Algebra Classes

// Distance algebra for forests
class DoubleDist_Algebra : public Algebra<double,RNA_Alphabet> {
private:
    Score s_;

public:
    double empty() const {
        return 0;
    };

    double replace(RNA_Alphabet a,double down, RNA_Alphabet b, double over) const {
        if (a==ALPHA_BASEPAIR && b == ALPHA_BASEPAIR)
            return s_.bp_rep_score_+down+over;
        else {
            if (a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
                return INT_MAX;
            else {
                if (a==b)
                    return s_.b_match_score_+down+over;
                else
                    return s_.b_rep_score_+down+over;
            }
        }
    };

    double del(RNA_Alphabet a,double down,double over) const {
        if (a==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double insert(double down,RNA_Alphabet b,double over) const {
        if (b==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double choice(double a, double  b) const {
        return std::min(a,b);
    };

    double worst_score() const {
        return INT_MAX;
    };

    DoubleDist_Algebra(const Score &s)
            : s_(s) {};
};


// Distance algebra for forests
class AffineDoubleDist_Algebra : public AlgebraAffine<double,RNA_Alphabet> {
private:
    Score s_;

public:
    double empty() const {
        return 0;
    };

    double replace(RNA_Alphabet a,double down, RNA_Alphabet b, double over) const {
        if (a==ALPHA_BASEPAIR && b == ALPHA_BASEPAIR)
            return s_.bp_rep_score_+down+over;
        else {
            if (a==ALPHA_BASEPAIR || b==ALPHA_BASEPAIR)
                return INT_MAX;
            else {
                if (a==b)
                    return s_.b_match_score_+down+over;
                else
                    return s_.b_rep_score_+down+over;
            }
        }
    };

    double del(RNA_Alphabet a,double down,double over) const {
        if (a==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double insert(double down,RNA_Alphabet b,double over) const {
        if (b==ALPHA_BASEPAIR)
            return s_.bp_indel_score_+down+over;
        else
            return s_.b_indel_score_+down+over;
    };

    double choice(double a, double  b) const {
        return std::min(a,b);
    };

    double worst_score() const {
        return INT_MAX;
    };

    double delO(RNA_Alphabet a,double down,double over) const {
        if (a==ALPHA_BASEPAIR)
            return s_.bp_indel_open_score_+down+over;
        else
            return s_.b_indel_open_score_+down+over;
    };

    double insertO(double down,RNA_Alphabet b,double over) const {
        if (b==ALPHA_BASEPAIR)
            return s_.bp_indel_open_score_+down+over;
        else
            return s_.b_indel_open_score_+down+over;
    };

    AffineDoubleDist_Algebra(const Score &s)
            : s_(s) {};
};


#endif
