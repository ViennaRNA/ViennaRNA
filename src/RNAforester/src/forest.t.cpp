#ifndef _FOREST_T_CPP_
#define _FOREST_T_CPP_

#include "misc.h"
#include "forest.h"
#include <cstring>
// Constructor and Destructor

template <class L>
Forest<L>::Forest(size_type size): ForestBase(size) {
    lb_=new L[size];
}

template <class L>
Forest<L>::Forest(const Forest<L> &f) : ForestBase(f) {
    lb_=new L[f.size()];

    // copy labels
    for (size_type i=0; i<size(); i++)
        lb_[i]=f.lb_[i];
}

template <class L>
Forest<L>::Forest(std::istream &s) {
    /*
    // read forest from file

    // first of all save the size
    s.read(static_cast<char*>(size_),sizeof(ForestBase::size_type));

    // initialize structures
    ::ForestBase(size_);
    lb_=new L[f.size()];

    // save the arrays
    s.read(static_cast<char*>(lb_),sizeof(label_type)*size_);
    s.read(static_cast<char*>(rb_),sizeof(ForestBase::size_type)*size_);
    s.read(static_cast<char*>(noc_),sizeof(ForestBase::size_type)*size_);
    s.read(static_cast<char*>(sumUpCSF_),sizeof(ForestBase::size_type)*size_);
    s.read(static_cast<char*>(rmb_),sizeof(ForestBase::size_type)*size_);
    */
}

template <class L>
Forest<L>::~Forest() {
    DELETE_ARRAY(lb_);
}

/* ****************************************** */
/*            Private functions               */
/* ****************************************** */

template<class L>
void Forest<L>::print(std::ostream &s,size_type node,bool brackets) const {
    if (brackets)
        s << "(";

    s <<  "{" << node << "}";
    showLabel(s,lb_[node]);
    if (noc_[node])
        print(s,node+1,true);
    if (rb_[node])
        print(s,rb_[node],false);

    if (brackets)
        s << ")";
}

// Calculate maximal score of csf (i,j) for a certain Algebra.
// It is assumed that a perfect match of the forest obtains the best
// possible score.

template <class L>
template <class R>
R Forest<L>::maxScore(const Algebra<R,L> &alg, size_type i, size_type j) const {
    R down, over;

    if (j==0)
        return 0;

    if (isLeaf(i)) {
        over=maxScore(alg,rb(i),j-1);
        return alg.replace(label(i),0,label(i),over);
    } else {
        down=maxScore(alg,i+1,noc(i));
        over=maxScore(alg,rb(i),j-1);
        return alg.replace(label(i),down,label(i),over);
    }
}

// Calculate maximal score of csf (i,j) for a certain RNA_Algebra.
// It is assumed that a perfect match of the forest obtains the best
// possible score.

template <class L>
template <class R>
R Forest<L>::maxScore(const RNA_Algebra<R,L> &alg, unsigned int i, unsigned int j) const {
    R down, over;

    if (j==0)
        return 0;

    if (isLeaf(i)) {
        over=maxScore(alg,rb(i),j-1);
        return alg.replace(label(i),0,label(i),over);
    } else {
        down=maxScore(alg,i+1+1,noc(i)-2);
        over=maxScore(alg,rb(i),j-1);
        return alg.replacepair(label(i+1),label(i+1),down,label(getRightmostBrotherIndex(i+1)),label(getRightmostBrotherIndex(i+1)),over);
    }
}

// Calculate maximal score of csf (i,j) for a certain Algebra.
// It is assumed that a perfect match of the forest obtains the best
// possible score.

template <class L>
template <class R>
R Forest<L>::maxScore(const AlgebraAffine<R,L> &alg, size_type i, size_type j) const {
    R down, over;

    if (j==0)
        return 0;

    if (isLeaf(i)) {
        over=maxScore(alg,rb(i),j-1);
        return alg.replace(label(i),0,label(i),over);
    } else {
        down=maxScore(alg,i+1,noc(i));
        over=maxScore(alg,rb(i),j-1);
        return alg.replace(label(i),down,label(i),over);
    }
}

// Calculate maximal score of csf (i,j) for a certain RNA_Algebra.
// It is assumed that a perfect match of the forest obtains the best
// possible score.

template <class L>
template <class R>
R Forest<L>::maxScore(const RNA_AlgebraAffine<R,L> &alg, unsigned int i, unsigned int j) const {
    R down, over;

    if (j==0)
        return 0;

    if (isLeaf(i)) {
        over=maxScore(alg,rb(i),j-1);
        return alg.replace(label(i),0,label(i),over);
    } else {
        down=maxScore(alg,i+1+1,noc(i)-2);
        over=maxScore(alg,rb(i),j-1);
        return alg.replacepair(label(i+1),label(i+1),down,label(getRightmostBrotherIndex(i+1)),label(getRightmostBrotherIndex(i+1)),over);
    }
}
/* ****************************************** */
/*            Protected function              */
/* ****************************************** */

template <class L>
void Forest<L>::initialize(size_type size, bool anchored) {
    ForestBase::initialize(size, anchored);
    lb_=new L[size];
		if (anchored) {
			anchors_=new size_type[size];
			memset(anchors_,0,sizeof(size_type)*size);
		}
}

/* ****************************************** */
/*             Public functions               */
/* ****************************************** */


template<class L>
void Forest<L>::printDot(std::ostream &s) const {
    size_type i,h;

    s << "digraph forest" << std::endl << "{" << std::endl;

    // edges
    for (i=0; i<size_; i++) {
        s << i;
        if (noc_[i]) {
            s << " -> {";
            h=i+1;
            for (unsigned int r=0; r<getMaxLength(i+1); r++) {
                s << h << " ";
                h=rb_[h];
            }
            s << "}";
        }
        s << std::endl;
    }

    // labels
    s << std::endl << std::endl;
    for (i=0; i<size_; i++) {
        s << i << "[label=\"";
        showLabel(s,lb_[i]);
        s << "\"]" << std::endl;
    }

    s << "}";
}

template<class L>
void Forest<L>::save(std::ostream &s) const {
    /*
    // save the pforest to stream in binary format

    // first of all save the size
    s.write(reinterpret_cast<char*>(&size_),sizeof(ForestBase::size_type));
    // save the arrays
    s.write(reinterpret_cast<char*>(lb_),sizeof(label_type)*size_);
    s.write(reinterpret_cast<char*>(rb_),sizeof(ForestBase::size_type)*size_);
    s.write(reinterpret_cast<char*>(noc_),sizeof(ForestBase::size_type)*size_);
    s.write(reinterpret_cast<char*>(sumUpCSF_),sizeof(ForestBase::size_type)*size_);
    s.write(reinterpret_cast<char*>(rmb_),sizeof(ForestBase::size_type)*size_);
    */
}



#endif








