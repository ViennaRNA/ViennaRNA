#ifndef _TREE_EDIT_T_CPP_
#define _TREE_EDIT_T_CPP_

#include <algorithm>
#include <cassert>

#include "alignment.h"
#include "debug.h"

#include "misc.t.cpp"
#include "forestsz.t.cpp"

template<class R,class L>
Mapping<R,L>::Mapping(const ForestSZ<L> *f1, const ForestSZ<L> *f2, const SZAlgebra<R,L> &alg) {
    assert(f1 != NULL);
    assert(f2 != NULL);

    unsigned long m,n,h,cols;
    unsigned long i,j,k,l,r,s;
    R score,h_score;


    // alloc space for the score matrices
    m=f1->size();
    n=f2->size();

    mtrxSizeTD_=(m)*(n);
    mtrxSizeFD_=(m+1)*(n+1);	// +1 because index 0 means the empty forest
    mtrxTD_=new R[mtrxSizeTD_];
    mtrxFD_=new R[mtrxSizeFD_];
    rowStartTD_=new unsigned long[m];
    rowStartFD_=new unsigned long[m+1];

    f1_ = new ForestSZ<L>(*f1);
    f2_ = new ForestSZ<L>(*f2);
    alg_=&alg;

    cols=n;
    rowStartTD_[0]=0;
    for (h=1; h<m; h++)
        rowStartTD_[h]=rowStartTD_[h-1]+cols;

    cols=n+1;
    rowStartFD_[0]=0;
    for (h=1; h<m+1; h++)
        rowStartFD_[h]=rowStartFD_[h-1]+cols;

    // calculate edit distance
    for (i=0; i<m; i++) {					// for all pairs of subtrees
        if (!f1_->keyroot(i))
            continue;

        for (j=0; j<n; j++) {
            if (!f2_->keyroot(j))
                continue;

            // calculate forest distance

            // edit to empty forest
            setMtrxFDVal(0,0,0);

            for (k=f1_->lml(i),r=1; k<=i; k++,r++)		// r is the indexpos of k
                setMtrxFDVal(r,0,alg_->del(f1_->label(k),getMtrxFDVal(r-1,0)));

            for (l=f2_->lml(j),s=1; l<=j; l++,s++)
                setMtrxFDVal(0,s,alg_->insert(getMtrxFDVal(0,s-1),f2_->label(l)));	// s is the indexpos of j

            for (k=f1_->lml(i),r=1; k<=i; k++,r++) {
                for (l=f2_->lml(j),s=1; l<=j; l++,s++) {
                    // fdist(k,i,l,j)
                    // lml(k)==lml(i) && lml(l)==lml(j)
                    if (f1_->lml(k)==f1_->lml(i) && f2_->lml(l)==f2_->lml(j)) {
                        score=alg_->replace(f1_->label(k),getMtrxFDVal(r-1,s-1),f2_->label(l));

                        h_score=alg_->del(f1_->label(k),getMtrxFDVal(r-1,s));
                        score=alg.choice(score,h_score);

                        h_score=alg_->insert(getMtrxFDVal(r,s-1),f2_->label(l));
                        score=alg.choice(score,h_score);

                        setMtrxFDVal(r,s,score);
                        setMtrxTDVal(k,l,score);
                    } else {
                        long p,q;
                        p=f1_->lml(k) - f1_->lml(i);
                        q=f2_->lml(l) - f2_->lml(j);

                        score=getMtrxFDVal(p,q) + getMtrxTDVal(k,l);

                        h_score=alg_->del(f1_->label(k),getMtrxFDVal(r-1,s));
                        score=alg.choice(score,h_score);

                        h_score=alg_->insert(getMtrxFDVal(r,s-1),f2_->label(l));
                        score=alg.choice(score,h_score);

                        setMtrxFDVal(r,s,score);
                    }
                }
            }
        }
    }

    // until here the original tree edit distance was calculated
    // to allow edit distances between forests, the distance for a virtual root node is calculated
    // leftmost leaf of the root node is the first node in postorder traversal
    // The distances are calculated until the last tree node in postorder traversal which then holds the result
    // for a virtual root node that is matched with zero cost

    // edit to empty forest
    setMtrxFDVal(0,0,0);

    for (k=0,r=1; k<m; k++,r++)		// r is the indexpos of k
        setMtrxFDVal(r,0,alg_->del(f1_->label(k),getMtrxFDVal(r-1,0)));

    for (l=0,s=1; l<n; l++,s++)
        setMtrxFDVal(0,s,alg_->insert(getMtrxFDVal(0,s-1),f2_->label(l)));	// s is the indexpos of j

    for (k=0,r=1; k<m; k++,r++) {
        for (l=0,s=1; l<n; l++,s++) {

            long p,q;
            p=f1_->lml(k);
            q=f2_->lml(l);

            score=getMtrxFDVal(p,q) + getMtrxTDVal(k,l);

            h_score=alg_->del(f1_->label(k),getMtrxFDVal(r-1,s));
            score=alg.choice(score,h_score);

            h_score=alg_->insert(getMtrxFDVal(r,s-1),f2_->label(l));
            score=alg.choice(score,h_score);

            setMtrxFDVal(r,s,score);
        }
    }

    optimum_=getMtrxFDVal(m,n);
}

template<class R,class L>
Mapping<R,L>::~Mapping() {
    DELETE(f1_);
    DELETE(f2_);

    delete[] mtrxTD_;
    delete[] mtrxFD_;
    delete[] rowStartTD_;
    delete[] rowStartFD_;
}
#endif
