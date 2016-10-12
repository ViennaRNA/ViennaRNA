#include "misc.h"
#include "forestbase.h"
#include "debug.h"

#include <algorithm>
#include <stdio.h>
#include <stddef.h>
#include <string.h>

#include "misc.t.cpp"

/* ****************************************** */
/*    Constructor and Destructor functions    */
/* ****************************************** */

// TODO selbst wissen ob anch
ForestBase::ForestBase(size_type size) {
    initialize(size,false);//not anchored
};

ForestBase::ForestBase(const ForestBase &fBase) {
    initialize(fBase.size(),false);//not anchored

    memcpy(rb_,fBase.rb_,sizeof(size_type)*size_);
    memcpy(noc_,fBase.noc_,sizeof(size_type)*size_);
    memcpy(sumUpCSF_,fBase.sumUpCSF_,sizeof(size_type)*(size_+1));
    memcpy(rmb_,fBase.rmb_,sizeof(size_type)*size_);

    isSumUpCSFValid_=fBase.isSumUpCSFValid_;
    isRMBValid_=fBase.isRMBValid_;
}

ForestBase::~ForestBase() {
	// TODO anch. deleten
    DELETE_ARRAY(rb_);
    DELETE_ARRAY(noc_);
    DELETE_ARRAY(sumUpCSF_);
    DELETE_ARRAY(rmb_);
};

/* ****************************************** */
/*            Private functions               */
/* ****************************************** */

ForestBase::size_type ForestBase::countLeaves(size_type i) const {
    size_type numLeaves=0;

    for (size_type k=0; k<i; k++) {
        if (isLeaf(k))
            numLeaves++;
    }

    return numLeaves;
};


/* ****************************************** */
/*            Protected function              */
/* ****************************************** */

void ForestBase::calcSumUpCSF() {
    // sum CSFs
    sumUpCSF_[0] = 1;
    for (size_type i=1; i<size_; i++) {
        sumUpCSF_[i] = sumUpCSF_[i-1]+getNumRightBrothers(i-1)+1;
    }

    sumUpCSF_[size_] = sumUpCSF_[size_-1]+1;  // last node is a leaf, so the number of siblings is one

    isSumUpCSFValid_=true;
}

void ForestBase::calcRMB() {
    for (int i=size_-1; i>=0; i--) {
        if (rb_[i])
            rmb_[i]=rmb_[rb_[i]];
        else
            rmb_[i]=i;
    }

    isRMBValid_=true;
}

void ForestBase::initialize(size_type size, bool anchored) {
    size_=size;
    rb_=new size_type[size];
    noc_=new size_type[size];
    sumUpCSF_=new size_type[size+1];
    rmb_=new size_type[size];
		if (anchored)
			anchors_ = new size_type[size];

    isSumUpCSFValid_=false;
    isRMBValid_=false;
}

/* ****************************************** */
/*             Public functions               */
/* ****************************************** */

ForestBase::size_type ForestBase::rb(size_type i, size_type k) const {
    if (k==0)
        return i;
    else
        return rb(rb(i),k-1);
}

ForestBase::size_type ForestBase::maxDegree() const {
    size_type degree=0;

    for (size_type i=0; i<size_; i++)
        degree=std::max(degree,getMaxLength(i));

    return degree;
}

ForestBase::size_type ForestBase::numLeaves() const {
    size_type count=0;

    for (size_type i=0; i<size_; i++)
        if (isLeaf(i))
            count++;

    return count;
}

ForestBase::size_type ForestBase::getNumRightBrothers(size_type node) const {
    size_type numRbs=0;

    while (rb_[node]) {
        numRbs++;
        node=rb_[node];
    }

    return numRbs;
}

ForestBase::size_type ForestBase::maxDepth() const {
    size_type *maxDepth;
    size_type m,h,r,i;

    maxDepth=new size_type[size_];

    for (r=0, i=size_-1; r<size_; r++,i--) {
        if (isLeaf(i))
            maxDepth[i]=1;
        else {
            m=0;
            h=i+1;		// h = first child
            do {
                m=std::max(m,maxDepth[h]);
                h=rb(h);
            } while (h);

            maxDepth[i]=1+m;
        }
    }

    // calculate maxDepth for the forest
    m=0;
    h=0;		// h = first root
    do {
        m=std::max(m,maxDepth[h]);
        h=rb(h);
    } while (h);

    delete maxDepth;
    return m;
}



