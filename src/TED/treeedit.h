#ifndef _TREE_EDIT_H_
#define _TREE_EDIT_H_

#include <cassert>
#include <vector>
#include "algebra.h"
#include "forestsz.h"

template<class R,class L>
class Mapping {
    R *mtrxTD_;
    R *mtrxFD_;
    unsigned long *rowStartTD_;
    unsigned long *rowStartFD_;
    unsigned long mtrxSizeTD_;
    unsigned long mtrxSizeFD_;
    R optimum_;

    ForestSZ<L> *f1_;
    ForestSZ<L> *f2_;
    const SZAlgebra<R,L> *alg_;

    inline R getMtrxTDVal(unsigned long i,unsigned long j) const {
        assert(rowStartTD_[i]+j<mtrxSizeTD_);
        return mtrxTD_[rowStartTD_[i]+j];
    };

    inline void setMtrxTDVal(unsigned long i,unsigned long j,R val) {
        assert(rowStartTD_[i]+j<mtrxSizeTD_);
        mtrxTD_[rowStartTD_[i]+j]=val;
    };

    inline R getMtrxFDVal(unsigned long i,unsigned long j) const {
        assert(rowStartFD_[i]+j<mtrxSizeFD_);
        return mtrxFD_[rowStartFD_[i]+j];
    };

    inline void setMtrxFDVal(unsigned long i,unsigned long j,R val) {
        assert(rowStartFD_[i]+j<mtrxSizeFD_);
        mtrxFD_[rowStartFD_[i]+j]=val;
    };

public:
    Mapping(const ForestSZ<L> *f1, const ForestSZ<L> *f2,const SZAlgebra<R,L> &alg);
    ~Mapping();

    R getGlobalOptimum() {
        return optimum_;
    };
};

#endif
