#ifndef _FOREST_H_
#define _FOREST_H_

#include <string>
#include <iostream>
#include <vector>

#include "algebra.h"
#include "misc.h"
#include "forestbase.h"



/** The template class Forest<L> is the basic class that is handled by the alignmnent
 *  algorithms of the RNAforester template library. The template parameter L is the
 *  datatype of the forest labels.
 */
template <class L>
class Forest : public ForestBase {
    template <class R, class LX, class AL>
    friend class Alignment;
		//friend std::ostream& operator<< <>(std::ostream &s, Forest<L> &f);
		friend std::ostream& operator<< (std::ostream &s, Forest<L> &f) {
			f.print(s, 0, true);
			return s;
		}

public:
    typedef typename ForestBase::size_type size_type;

private:
    /** internal function used by '<<' operator to print bracket notation of the forest*/
    void print(std::ostream &s,size_type node,bool brackets) const;

    /** Calculate maximal score of csf (i,j) for a certain Algebra.
    *  It is assumed that a perfect match of the forest obtains the best
    *  possible score.
    */
    template <class R>
    R maxScore(const AlgebraAffine<R,L> &alg, size_type i, size_type j) const;
    template <class R>
    R maxScore(const RNA_AlgebraAffine<R,L> &alg, size_type i, size_type j) const;
    template <class R>
    R maxScore(const Algebra<R,L> &alg, size_type i, size_type j) const;
    template <class R>
    R maxScore(const RNA_Algebra<R,L> &alg, size_type i, size_type j) const;

    /** Function showLabel is used by print routines of Forest */
//	virtual void showLabel(std::ostream &s,char c) const {s << c;};
    virtual void showLabel(std::ostream &s,L l) const {
        s << 'X';
    };

    /** Function makeLabel is used by function buildForest */
//	virtual void makeLabel(char &a,char c) {a=c;};

    void makeLabel(L &a,char c) {};


protected:
    L *lb_;


public:
    typedef L label_type;

    /** Default Constructor */
    Forest() : ForestBase(), lb_(NULL) {};
    /** Construct a Forest with 'size' nodes */
    Forest(size_type size);
    /** Copy constructor */
    Forest(const Forest<L> &f);
    /** Read Forest from binary file */
    Forest(std::istream &s);

    virtual ~Forest();

    /** returns label of node i */
    inline L label(size_type i) const {
        return lb_[i];
    };

    /** Calculate maximal score of a forest alignment against itself for a certain Algebra.
    *  It is assumed that a perfect match of the forest obtains the best
    *  possible score.
    */
    template <class R>
    inline R maxScore(const Algebra<R,L> &alg) const {return maxScore(alg,0,getMaxLength(0));}

    template <class R>
    inline R maxScore(const RNA_Algebra<R,L> &alg) const { return maxScore(alg,0,getMaxLength(0)); }

    template <class R>
    inline R maxScore(const AlgebraAffine<R,L> &alg) const {return maxScore(alg,0,getMaxLength(0));}

    template <class R>
    inline R maxScore(const RNA_AlgebraAffine<R,L> &alg) const { return maxScore(alg,0,getMaxLength(0)); }

		// TODO war nicht immer public, scheint aber irgendwie noetig zu sein
    void initialize(size_type size, bool anchored);

    /** Print forest in GraphViz format */
    void printDot(std::ostream &s) const;

    /** Save forest to binary file */
    void save(std::ostream &s) const;

void printMembers() const {
	std::cout << std::endl;
	std::cout << "size " << size_ << std::endl;

	std::cout<<"rb"<<std::endl;
	for (size_type i=0;i<size();i++) {
		std::cout << rb_[i] << ",";
	}

	std::cout << "noc" << std::endl;
	for (size_type i=0;i<size();i++) {
		std::cout << noc_[i] << ",";
	}
	std::cout << std::endl;

	std::cout << "sumUpCSF" << std::endl;
	for (size_type i=0;i<size();i++) {
		std::cout << sumUpCSF_[i] << ",";
	}
	std::cout << std::endl;

	std::cout << "rmb" << std::endl;
	for (unsigned int i = 0; i <= size_-1; i++) {
		std::cout << rmb_[i] << ",";
	}
	std::cout << std::endl;

	std::cout<<"anchors"<<std::endl;
	for (unsigned int i = 0; i <= size_-1; i++) {
;//		std::cout << anchors_[i] << ",";
	}
	std::cout << std::endl;

	std::cout<<"lb"<<std::endl;
	for (size_type i=0;i<size();i++) {
		//std::cout << lb_[i] << ",";
		showLabel(std::cout, lb_[i]);
	}
	std::cout<<std::endl;

}



};

/** Stream output operator for Forest<L> */
template <class L>
std::ostream& operator<< (std::ostream &s, Forest<L> &f);

#endif

