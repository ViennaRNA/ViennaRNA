#ifndef _ALIGNMENT_T_CPP_
#define _ALIGNMENT_T_CPP_

#include "alignment.h"
#include "forest.t.cpp"
#include <fstream>

// Constructor and Destructor

// TODO merge constructors!
// TODO divide into fill und backtrace functionality
// TODO pIndels anschliessen

// costructor for general and profile alignments
template<class R,class L,class AL>
Alignment<R,L,AL>::Alignment(const Forest<L> *f1, const Forest<L> *f2, const bool topdown, const bool anchored, const bool printBacktrace)
        : f1_(f1),
        f2_(f2),
				topdown_(topdown),// kann weg als memeber nicht noetig
				anchored_(anchored),
				pIndels_(true),
				printBacktrace_(printBacktrace),
        suboptimalsPercent_(100) {

}

template<class R,class L,class AL>
Alignment<R,L,AL>::~Alignment() {
    //TODO in unterklassen? delete mtrx_;
}


// CYK style parser / bottom up matrix fill 
template<class R, class L, class AL>
void Alignment<R, L, AL>::fillMatricesBottomUp(bool local, bool RNA, bool speedup) {
      makeFirstCell();
      makeFirstRow();
      makeFirstCol();
      if (local)
        makeInnersLocal(RNA, speedup);
      else
        makeInnersGlobal(RNA, speedup);// hier wird weniger ausgerechnet
}

// Unger style parser / bottom up matrix fill 
template<class R, class L, class AL>
void Alignment<R, L, AL>::fillMatricesTopDown(bool local, bool RNA, bool speedup) {

  // recursion start values
  unsigned int i = 0, j = this->f1_->getMaxLength(0);
  unsigned int k = 0, l = this->f2_->getMaxLength(0);
  CSFPair startpair = CSFPair(i,j,k,l);
  //this->calls = 0;

  if (anchored_) {
    recursiveFillAnchored(startpair, RNA, speedup);
  }
  else {
    recursiveFill(startpair, RNA, speedup);
  }
  //std::cout << "all = " << this->f1_.getNumCSFs()*this->f2_.getNumCSFs() << std::endl;
  //std::cout << "calls = " << calls << std::endl;
  //std::cout << "perc = " << (((double) calls)/((double) this->f1_.getNumCSFs()*this->f2_.getNumCSFs())) << std::endl;
}

template<class R, class L, class AL>
void Alignment<R, L, AL>::recursiveFill(CSFPair p, bool RNA, bool speedup) {

  unsigned int r = 0, h = 0;
  //this->calls++;

  // easiest: already done;
  if (computed(this->f1_->indexpos(p.i,p.j), this->f2_->indexpos(p.k, p.l)))
    return;

  // no recursion for empty ali
  // -> jump straight to computation

  // replacement
  if (p.j > 0 && p.l > 0) {
		if (RNA && this->f1_->down(p.i) && this->f2_->down(p.k)) {
			if (! computed(this->f1_->mdown(p.i),this->f2_->mdown(p.k)))
				recursiveFill(CSFPair(p.i + 1 + 1, this->f1_->noc(p.i) - 2, p.k + 1 + 1, this->f2_->noc(p.k) - 2), RNA, speedup);
		}
		else {
			if (! computed(this->f1_->down(p.i), this->f2_->down(p.k)))
				recursiveFill(CSFPair(p.i + 1, this->f1_->noc(p.i), p.k + 1, this->f2_->noc(p.k)), RNA, speedup);
		}
		if (! computed(this->f1_->over(p.i, p.j) , this->f2_->over(p.k, p.l)))
			recursiveFill(CSFPair(this->f1_->rb(p.i), p.j - 1, this->f2_->rb(p.k), p.l - 1), RNA, speedup);
  }
  // deletion
  if (p.j > 0) {
    if (this->f1_->noc(p.i)==0 && speedup) { // no child
			if (! computed(this->f1_->over(p.i,p.j),this->f2_->indexpos(p.k,p.l))) {
				recursiveFill(CSFPair(this->f1_->rb(p.i), p.j - 1, p.k, p.l), RNA, speedup);
			}
    }
    else  if (this->f1_->rb(p.i)==0 && speedup) { // no right brother
			if (! computed(this->f1_->down(p.i),this->f2_->indexpos(p.k,p.l))){
				recursiveFill(CSFPair(p.i + 1, this->f1_->noc(p.i), p.k, p.l), RNA, speedup);
			}
    }
		else {
			h = p.k; // h is the node where the suffix of the split begins
			for (r = 0; r <= p.l; r++) {// for all splits of f2
				if (! computed(this->f1_->down(p.i) , this->f2_->indexpos(p.k, r)))
					recursiveFill(CSFPair(p.i + 1, this->f1_->noc(p.i), p.k, r), RNA, speedup);
				if (! computed(this->f1_->over(p.i, p.j) , this->f2_->indexpos(h, p.l - r)))
					recursiveFill(CSFPair(this->f1_->rb(p.i), p.j - 1, h, p.l - r), RNA, speedup);
				h = this->f2_->rb(h);
			}
		}
  }
   // insert
  if (p.l > 0) {
    if (this->f2_->noc(p.k)==0 && speedup) { // no child
			if (! computed(this->f1_->indexpos(p.i,p.j),this->f2_->over(p.k,p.l))) {
				recursiveFill(CSFPair(p.i, p.j,this->f2_->rb(p.k), p.l - 1), RNA, speedup);
			}
    }
    else if (this->f2_->rb(p.k)==0 && speedup) { // no right brother
			if (! computed(this->f1_->indexpos(p.i,p.j),this->f2_->down(p.k))) {
				recursiveFill(CSFPair(p.i, p.j, p.k + 1, this->f2_->noc(p.k)), RNA, speedup);
			}
    }
		else {
			h = p.i;
			for (r = 0; r <= p.j; r++) {// for all splits of f1
				if (! computed(this->f1_->indexpos(p.i, r) , this->f2_->down(p.k)))
					recursiveFill(CSFPair(p.i, r, p.k + 1, this->f2_->noc(p.k)), RNA, speedup);
				if (! computed(this->f1_->indexpos(h, p.j - r) , this->f2_->over(p.k, p.l)))
					recursiveFill(CSFPair(h, p.j - r, this->f2_->rb(p.k), p.l - 1), RNA, speedup);
				h = this->f1_->rb(h);
			}
		}
  }

  // finally: computation of the entry
  compareCSFPairTD(p, speedup);
	return;
}


template<class R,class L,class AL>
AlignmentLinear<R,L,AL>::AlignmentLinear(const Forest<L> *f1, const Forest<L> *f2, const Algebra<R,L> &alg, const bool topdown, const bool anchored, const bool local, const bool printBacktrace, bool speedup)
        : Alignment<R,L,AL>(f1,f2,topdown,anchored,printBacktrace) {

    // alloc space for the score matrix, backtrace structure,
    // and , if wanted, for the calculation-order-matrix
		mtrx_ = new TAD_DP_TableLinear<R>(this->f1_->getNumCSFs(),this->f2_->getNumCSFs(),alg.worst_score());
    // initialize variables
    alg_ = &alg;
    rnaAlg_ = NULL;
    this->localOptimum_ = alg.worst_score();
    // align forests f1_ and f2_
		bool RNA = false;
		if (topdown)
			this->fillMatricesTopDown(local, RNA, speedup);
		else
			this->fillMatricesBottomUp(local, RNA, speedup);
}


// constructor for RNA alignments
template<class R,class L,class AL>
AlignmentLinear<R,L,AL>::AlignmentLinear(const Forest<L> *f1, const Forest<L> *f2, const RNA_Algebra<R,L> &rnaAlg, const bool topdown, const bool anchored, const bool local, const bool printBacktrace, bool speedup)
        : Alignment<R,L,AL>(f1,f2,topdown,anchored,printBacktrace) {

    // alloc space for the score matrix, backtrace structure and,
    // if wanted, for the calculation-order-matrix
		mtrx_ = new TAD_DP_TableLinear<R>(this->f1_->getNumCSFs(),this->f2_->getNumCSFs(),rnaAlg.worst_score());
    // initialize variables
    rnaAlg_ = &rnaAlg;
    alg_ = (const Algebra<R,L>*)&rnaAlg;
    this->localOptimum_ = rnaAlg.worst_score();
    // align forests f1 and f2
		bool RNA = true;
		if (topdown)
			this->fillMatricesTopDown(local, RNA, speedup);
		else
			this->fillMatricesBottomUp(local, RNA, speedup);
}

// Private functions


template<class R, class L, class AL>
void AlignmentLinear<R, L, AL>::makeFirstCell() {
    setMtrxVal(0, 0, alg_->empty());
}

template<class R, class L, class AL>
void AlignmentLinear<R, L, AL>::makeFirstRow() {
    // align f2 to the empty forest (fill first row of array)
    // for all nodes in f1
    unsigned long n = this->f2_->size();
    for (long k = n - 1; k >= 0; k--) {
        // for all non empty csfs induced by k
        for (unsigned int l = 1; l <= this->f2_->getMaxLength(k); l++) {
            R score = alg_->insert(getMtrxVal(0, this->f2_->down(k)),
                                 this->f2_->label(k),
                                 getMtrxVal(0, this->f2_->over(k, l)));
            setMtrxVal(0, this->f2_->indexpos(k, l), score);
        }
    }
}

template<class R, class L, class AL>
void AlignmentLinear<R, L, AL>::makeFirstCol() {
    // align f1 to the empty forest (fill first column of array)
    // for all nodes in f1
    unsigned long m = this->f1_->size();
    for (long i = m - 1; i >= 0; i--) {
        // for all non empty csfs induced by i
        for (unsigned int j = 1; j <= this->f1_->getMaxLength(i); j++) {
            R score = alg_->del(this->f1_->label(i),
                              getMtrxVal(this->f1_->down(i), 0),
                              getMtrxVal(this->f1_->over(i, j), 0));
            setMtrxVal(this->f1_->indexpos(i, j), 0, score);
        }
    }
}

template<class R,class L,class AL>
void Alignment<R,L,AL>::makeInnersLocal(bool RNA, bool speedup) {

    // align the rest
    unsigned long m = this->f1_->size();
    unsigned long n = this->f2_->size();
    for (long i=m-1; i>=0; i--) // for all nodes in f1_
        for (long k=n-1; k>=0; k--) // for all nodes in f1_
            for (unsigned int j=1; j<=this->f1_->getMaxLength(i); j++) // for all non empty csfs induced by i
                for (unsigned int l=1; l<=this->f2_->getMaxLength(k); l++) { // for all non empty csfs induced by k
                    compareCSFPair(CSFPair(i, j, k, l), RNA, speedup);
								}

    resetOptLocalAlignment(100);
}

template<class R,class L,class AL>
void Alignment<R,L,AL>::makeInnersGlobal(bool RNA, bool speedup) {

    // align the rest
    unsigned long m = this->f1_->size();
    unsigned long n = this->f2_->size();
    for (long i=m-1; i>=0; i--)  // for all nodes in f1_
        for (long k=n-1; k>=0; k--) { // for all nodes in f1_
            // compute matrix cols for subforests of global ali
            unsigned int l = this->f2_->getMaxLength(k);
            for (unsigned int j=1; j<=this->f1_->getMaxLength(i); j++) { // for all non empty csfs induced by i
                compareCSFPair(CSFPair(i,j,k,l), RNA, speedup);
                if (RNA && l>1)
                    compareCSFPair(CSFPair(i,j,k,l-1), RNA, speedup);
            }
            // compute matrix rows for subforests of global ali
            unsigned int j = this->f1_->getMaxLength(i);
            for (unsigned int l=1; l<=this->f2_->getMaxLength(k); l++) { // for all non empty csfs induced by k
                compareCSFPair(CSFPair(i,j,k,l), RNA, speedup);
                if (RNA && j>1)
                    compareCSFPair(CSFPair(i,j-1,k,l), RNA, speedup);
            }
        }
    resetOptLocalAlignment(100);
}

template<class R, class L, class AL>
void AlignmentLinear<R, L, AL>::compareCSFPairTD(CSFPair p, bool speedup) {
	//std::cout << "topdown" << std::endl;
  R score = alg_->worst_score();
  R h_score = alg_->worst_score();

  // empty alignment
  if (p.j == 0 && p.l == 0) {
    makeFirstCell();
  	setComputed(this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l));
    return;
  }
  // replacement
  if (p.j > 0 && p.l > 0) {
    // basepair replacement
    if (rnaAlg_ && this->f1_->down(p.i) && this->f2_->down(p.k)) {
        // must be two P nodes !!
        score = rnaAlg_->replacepair(this->f1_->label(p.i+1),
                                       this->f2_->label(p.k+1),
                                       getMtrxVal(this->f1_->mdown(p.i),this->f2_->mdown(p.k)),
                                       this->f1_->label(this->f1_->getRightmostBrotherIndex(p.i+1)),
                                       this->f2_->label(this->f2_->getRightmostBrotherIndex(p.k+1)),
                                       getMtrxVal(this->f1_->over(p.i,p.j),this->f2_->over(p.k,p.l)));
    } // replacement
    else {
        score = alg_->replace(this->f1_->label(p.i),
                                   getMtrxVal(this->f1_->down(p.i),this->f2_->down(p.k)),
                                   this->f2_->label(p.k),
                                   getMtrxVal(this->f1_->over(p.i,p.j),this->f2_->over(p.k,p.l)));
    }
  }
  // delete
  if (p.j > 0) {
    if (this->f1_->noc(p.i)==0 && speedup) { // no child
        h_score = alg_->del(this->f1_->label(p.i),0, getMtrxVal(this->f1_->over(p.i,p.j),this->f2_->indexpos(p.k,p.l)));
        score = alg_->choice(score,h_score);
    }
    else  if (this->f1_->rb(p.i)==0 && speedup) { // no right brother
        h_score = alg_->del(this->f1_->label(p.i), getMtrxVal(this->f1_->down(p.i),this->f2_->indexpos(p.k,p.l)), 0);
        score = alg_->choice(score,h_score);
    }
    else {
        unsigned long h = p.k;               // h is the node where the suffix of the split begins
        for (unsigned int r=0; r<=p.l; r++) { // for all splits of f2_
            h_score = alg_->del(this->f1_->label(p.i),
                                getMtrxVal(this->f1_->down(p.i),this->f2_->indexpos(p.k,r)),
                                getMtrxVal(this->f1_->over(p.i,p.j),this->f2_->indexpos(h,p.l-r)));
            score = alg_->choice(score,h_score);
            h = this->f2_->rb(h);
        }
    }
  }
  // insert
  if (p.l > 0) {
    if (this->f2_->noc(p.k)==0 && speedup) { // no child
        h_score = alg_->insert(0, this->f2_->label(p.k), getMtrxVal(this->f1_->indexpos(p.i,p.j),this->f2_->over(p.k,p.l)));
        score = alg_->choice(score,h_score);
    }
    else if (this->f2_->rb(p.k)==0 && speedup) { // no right brother
        h_score = alg_->insert(getMtrxVal(this->f1_->indexpos(p.i,p.j),this->f2_->down(p.k)), this->f2_->label(p.k), 0);
        score = alg_->choice(score,h_score);
    }
    else {
        unsigned long h = p.i;
        for (unsigned int r=0; r<=p.j; r++) { // for all splits of f1_
            h_score = alg_->insert(getMtrxVal(this->f1_->indexpos(p.i,r),this->f2_->down(p.k)),
                                   this->f2_->label(p.k),
                                   getMtrxVal(this->f1_->indexpos(h,p.j-r),this->f2_->over(p.k,p.l)));
            score = alg_->choice(score,h_score);
            h = this->f1_->rb(h);
        }
    }
  }
	// set value
  setMtrxVal(this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l), score);
  setComputed(this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l));
}

template<class R, class L, class AL>
void AlignmentLinear<R, L, AL>::compareCSFPair(CSFPair p, bool RNA, bool speedup) {
		R score;
    R h_score;
    // basepair replacement
    if (rnaAlg_ && this->f1_->down(p.i) && this->f2_->down(p.k)) {
        // must be two P nodes !!
        score = rnaAlg_->replacepair(this->f1_->label(p.i+1),
                                       this->f2_->label(p.k+1),
                                       getMtrxVal(this->f1_->mdown(p.i),this->f2_->mdown(p.k)),
                                       this->f1_->label(this->f1_->getRightmostBrotherIndex(p.i+1)),
                                       this->f2_->label(this->f2_->getRightmostBrotherIndex(p.k+1)),
                                       getMtrxVal(this->f1_->over(p.i,p.j),this->f2_->over(p.k,p.l)));
    } // replacement
    else {
        score = alg_->replace(this->f1_->label(p.i),
                                   getMtrxVal(this->f1_->down(p.i),this->f2_->down(p.k)),
                                   this->f2_->label(p.k),
                                   getMtrxVal(this->f1_->over(p.i,p.j),this->f2_->over(p.k,p.l)));
    }
    // delete
    if (this->f1_->noc(p.i)==0 && speedup) { // no child
        h_score = alg_->del(this->f1_->label(p.i),0, getMtrxVal(this->f1_->over(p.i,p.j),this->f2_->indexpos(p.k,p.l)));
        score = alg_->choice(score,h_score);
    } 
		else  if (this->f1_->rb(p.i)==0 && speedup) { // no right brother
        h_score = alg_->del(this->f1_->label(p.i), getMtrxVal(this->f1_->down(p.i),this->f2_->indexpos(p.k,p.l)), 0);
        score = alg_->choice(score,h_score);
    } 
		else {
        unsigned long h = p.k;               // h is the node where the suffix of the split begins
        for (unsigned int r=0; r<=p.l; r++) { // for all splits of f2_
            h_score = alg_->del(this->f1_->label(p.i),
                                getMtrxVal(this->f1_->down(p.i),this->f2_->indexpos(p.k,r)),
                                getMtrxVal(this->f1_->over(p.i,p.j),this->f2_->indexpos(h,p.l-r)));
            score = alg_->choice(score,h_score);
            h = this->f2_->rb(h);
        }
    }
    // insert
    if (this->f2_->noc(p.k)==0 && speedup) { // no child
        h_score = alg_->insert(0, this->f2_->label(p.k), getMtrxVal(this->f1_->indexpos(p.i,p.j),this->f2_->over(p.k,p.l)));
        score = alg_->choice(score,h_score);
    } 
		else if (this->f2_->rb(p.k)==0 && speedup) { // no right brother
        h_score = alg_->insert(getMtrxVal(this->f1_->indexpos(p.i,p.j),this->f2_->down(p.k)), this->f2_->label(p.k), 0);
        score = alg_->choice(score,h_score);
    } 
		else {
        unsigned long h = p.i;
        for (unsigned int r=0; r<=p.j; r++) { // for all splits of f1_
            h_score = alg_->insert(getMtrxVal(this->f1_->indexpos(p.i,r),this->f2_->down(p.k)),
                                   this->f2_->label(p.k),
                                   getMtrxVal(this->f1_->indexpos(h,p.j-r),this->f2_->over(p.k,p.l)));
            score = alg_->choice(score,h_score);
            h = this->f1_->rb(h);
        }
    }
    // set value
    setMtrxVal(this->f1_->indexpos(p.i,p.j),this->f2_->indexpos(p.k,p.l),score);
}

template<class R,class L,class AL>
void AlignmentLinear<R,L,AL>::recursiveFillAnchored(CSFPair p, bool RNA, bool speedup) {
    R score = this->alg_->worst_score();
    R h_score = this->alg_->worst_score();
    int a = this->f1_->getAnchor(p.i);
    int b = this->f2_->getAnchor(p.k);
    //this->calls++;

    // easiest: already done;
    if (computed(this->f1_->indexpos(p.i,p.j) , this->f2_->indexpos(p.k, p.l)))
        return;

    // empty alignment
    // no recursion for empty ali
    if (p.j == 0 && p.l == 0) {
        this->makeFirstCell();
        return;
    }

    // replacement
    // a and b both anchors or both nonanchors
    if (a==b && p.j > 0 && p.l > 0) {
        if (a>0 && b>0) {
            //std::cout << "Matched anchor " << a << "." << std::endl;
        }
    		if (RNA && this->f1_->down(p.i) && this->f2_->down(p.k)) {
						if (! computed(this->f1_->mdown(p.i),this->f2_->mdown(p.k)))
							recursiveFillAnchored(CSFPair(p.i + 1 + 1, this->f1_->noc(p.i) - 2, p.k + 1 + 1, this->f2_->noc(p.k) - 2), RNA, speedup);
        }
				else {
						if (! computed(this->f1_->down(p.i), this->f2_->down(p.k)))
							recursiveFillAnchored(CSFPair(p.i + 1, this->f1_->noc(p.i), p.k + 1, this->f2_->noc(p.k)), RNA, speedup);
				}
				if (! computed(this->f1_->over(p.i, p.j) , this->f2_->over(p.k, p.l)))
					recursiveFillAnchored(CSFPair(this->f1_->rb(p.i), p.j - 1, this->f2_->rb(p.k), p.l - 1), RNA, speedup);

				R ovr = getMtrxVal(this->f1_->over(p.i,p.j),this->f2_->over(p.k,p.l));
				// basepair replacement
				// must be two P nodes !!
				if (rnaAlg_ && this->f1_->down(p.i) && this->f2_->down(p.k)) {
					R dwn = getMtrxVal(this->f1_->mdown(p.i), this->f2_->mdown(p.k));
			    score = rnaAlg_->replacepair(this->f1_->label(p.i+1), this->f2_->label(p.k+1), dwn,
			                 this->f1_->label(this->f1_->getRightmostBrotherIndex(p.i+1)),
			                 this->f2_->label(this->f2_->getRightmostBrotherIndex(p.k+1)), ovr);
				}// base replacement
				else {
					R dwn = getMtrxVal(this->f1_->down(p.i),this->f2_->down(p.k));
			    score = alg_->replace(this->f1_->label(p.i), dwn, this->f2_->label(p.k), ovr);
			  }

				// replacement ende
				// two anchors were matched - done
        if (a>0 && b>0) {
          // set value
          setMtrxVal(this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l), score);
          setComputed(this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l));
          return;
        }
		}
		
    // deletion
    if (a==0 && p.j > 0) {
	
					L lbl = this->f1_->label(p.i);
			    if (this->f1_->noc(p.i)==0 && speedup) { // no child
							// right alignment
							if (! computed(this->f1_->over(p.i, p.j) , this->f2_->indexpos(p.k, p.l)))
								recursiveFillAnchored(CSFPair(this->f1_->rb(p.i), p.j - 1, p.k, p.l), RNA, speedup);
							R ovr = getMtrxVal(this->f1_->over(p.i,p.j),this->f2_->indexpos(p.k,p.l));
							h_score = alg_->del(lbl, 0, ovr);
			        score = alg_->choice(score,h_score);
			    } 
					else  if (this->f1_->rb(p.i)==0 && speedup) { // no right brother
							// down alignment
							if (! computed(this->f1_->down(p.i) , this->f2_->indexpos(p.k, p.l)))
								recursiveFillAnchored(CSFPair(p.i + 1, this->f1_->noc(p.i), p.k, p.l), RNA, speedup);
							R dwn = getMtrxVal(this->f1_->down(p.i),this->f2_->indexpos(p.k,p.l));
							h_score = alg_->del(lbl, dwn, 0);
							score = alg_->choice(score,h_score);
			    } 
					else {
						// h = node where suffix of split begins
						unsigned int h = p.k;
						// for all splits of f2
				  	for (unsigned int r = 0; r <= p.l; r++) {
								// down alignment
								if (! computed(this->f1_->down(p.i) , this->f2_->indexpos(p.k, r)))
									recursiveFillAnchored(CSFPair(p.i + 1, this->f1_->noc(p.i), p.k, r), RNA, speedup);
								// right alignment
								if (! computed(this->f1_->over(p.i, p.j) , this->f2_->indexpos(h, p.l - r)))
									recursiveFillAnchored(CSFPair(this->f1_->rb(p.i), p.j - 1, h, p.l - r), RNA, speedup);
								R dwn = getMtrxVal(this->f1_->down(p.i), this->f2_->indexpos(p.k, r));
								R ovr = getMtrxVal(this->f1_->over(p.i, p.j), this->f2_->indexpos(h, p.l - r));
						    h_score = alg_->del(lbl, dwn, ovr);
						      
				  	    score = alg_->choice(score, h_score);
								h = this->f2_->rb(h);
						}
					}
    }

    // insertion
    if (b==0 && p.l > 0) {
	
					L lbl = this->f2_->label(p.k);
			    if (this->f2_->noc(p.k)==0 && speedup) { // no child
							// right alignment
							if (! computed(this->f1_->indexpos(p.i, p.j) , this->f2_->over(p.k, p.l)))
								recursiveFillAnchored(CSFPair(p.i, p.j, this->f2_->rb(p.k), p.l - 1), RNA, speedup);
							R ovr = getMtrxVal(this->f1_->indexpos(p.i,p.j),this->f2_->over(p.k,p.l));
							h_score = alg_->insert(0, lbl, ovr);
			        score = alg_->choice(score,h_score);
			    } 
					else if (this->f2_->rb(p.k)==0 && speedup) { // no right brother
							// down alignment
							if (! computed(this->f1_->indexpos(p.i, p.j) , this->f2_->down(p.k)))
								recursiveFillAnchored(CSFPair(p.i, p.j, p.k + 1, this->f2_->noc(p.k)), RNA, speedup);
							R dwn = getMtrxVal(this->f1_->indexpos(p.i,p.j),this->f2_->down(p.k));
							h_score = alg_->insert(dwn, lbl, 0);
			        score = alg_->choice(score,h_score);
			    } 
					else {
							// h = node where suffix of split begins
						  unsigned int h = p.i;
						  // for all splits of f1
						  for (unsigned int r = 0; r <= p.j; r++) {
									// down alignment
									if (! computed(this->f1_->indexpos(p.i, r) , this->f2_->down(p.k)))
										recursiveFillAnchored(CSFPair(p.i, r, p.k + 1, this->f2_->noc(p.k)), RNA, speedup);
									// right alignment
									if (! computed(this->f1_->indexpos(h, p.j - r) , this->f2_->over(p.k, p.l)))
										recursiveFillAnchored(CSFPair(h, p.j - r, this->f2_->rb(p.k), p.l - 1), RNA, speedup);
									R dwn = getMtrxVal(this->f1_->indexpos(p.i, r), this->f2_->down(p.k));
									R ovr = getMtrxVal(this->f1_->indexpos(h, p.j - r), this->f2_->over(p.k, p.l));
						      h_score =  alg_->insert(dwn, lbl, ovr);
									//std::cout << "dwn " << dwn << " ovr " << ovr << " open " << open << std::endl;
									//std::cout << "score t " << score[t] << " h_score t " << h_score[t] << std::endl;
						      score = alg_->choice(score, h_score);
									h = this->f1_->rb(h);
						  }
			    }
    }
		if (score > 147483580) {
			std::cout << "score overflow problem - seen score " << std::endl;
			std::cout << score << " in matrix " << std::endl;
			exit(0);
		}
    // set value
    setMtrxVal(this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l), score);
    setComputed(this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l));
};

template<class R,class L,class AL>
unsigned int AlignmentLinear<R,L,AL>::backtrace(ForestAli<L,AL> &f, CSFPair p, unsigned int &node, int t) {

		if (this->printBacktrace_)
			std::cout << "backtracing from i=" << p.i << ", j=" << p.j << ", k=" << p.k << ", l=" << p.l << std::endl;
    // empty alignment
    if (p.j==0 && p.l==0) {
        return 0;
    }

    unsigned int p_node,rb_node,noc,rbs,r,h;
    R score = getMtrxVal(this->f1_->indexpos(p.i,p.j),this->f2_->indexpos(p.k,p.l));
    R h_score;
    p_node = node;
    node++;

		//std::cout << "score " << score << std::endl;
    // could it be a replacement
    if (p.j>0 && p.l>0) {
        // check for basepair replacement only if Algebra is of type RNA_Algebra
        if (rnaAlg_ && this->f1_->down(p.i) && this->f2_->down(p.k)) {
            h_score = rnaAlg_->replacepair(this->f1_->label(p.i+1), this->f2_->label(p.k+1),
                      getMtrxVal(this->f1_->mdown(p.i),this->f2_->mdown(p.k)), this->f1_->label(this->f1_->getRightmostBrotherIndex2(p.i+1)),
                      this->f2_->label(this->f2_->getRightmostBrotherIndex2(p.k+1)),
                      getMtrxVal(this->f1_->over(p.i,p.j),this->f2_->over(p.k,p.l)));

            if (score == h_score) {
								int anchor = 0;
								if (this->anchored_) {
									if (this->f1_->getAnchor(p.i)==this->f2_->getAnchor(p.k)) {
										anchor = this->f1_->getAnchor(p.i);
									}
									else {
										std::cout << "apparently matchend unequal anchors." << std::endl;
										exit(EXIT_FAILURE);
									}
								}
                // it is a basepair replacement
                f.makeRepLabel(p_node,this->f1_->label(p.i),this->f2_->label(p.k),anchor);   // P
                // set labels of leftmost child
                f.makeRepLabel(p_node+1,this->f1_->label(p.i+1),this->f2_->label(p.k+1),0);
                f.setRightBrotherIndex(p_node+1,p_node+1+1);
                f.setNumChildren(p_node+1,0);  // base node has no children
                node++;
                // down alignment
                assert(this->f1_->noc(p.i)>=2);
                assert(this->f2_->noc(p.k)>=2);
                noc = backtrace(f, CSFPair(p.i+1+1,this->f1_->noc(p.i)-2,p.k+1+1,this->f2_->noc(p.k)-2), node);  // !! mdown !!
                f.setNumChildren(p_node,noc+2);

                if (noc==0) {
                    f.setRightBrotherIndex(p_node+1,p_node+1+1);
                    f.setRightBrotherIndex(p_node+1+1,0);
                } else {
                    f.setRightBrotherIndex(f.getRightmostBrotherIndex2(p_node+1+1),node);
                }
                // set labels of leftmost child
                f.makeRepLabel(node,this->f1_->label(this->f1_->getRightmostBrotherIndex2(p.i+1+1)),this->f2_->label(this->f2_->getRightmostBrotherIndex2(p.k+1+1)),0);
                f.setRightBrotherIndex(node,0);
                f.setNumChildren(node,0);  // base node has no children
                node++;
                rb_node = node; // !!
                // right alignment
                rbs = backtrace(f, CSFPair(this->f1_->rb(p.i),p.j-1,this->f2_->rb(p.k),p.l-1), node);
                if (rbs)
                    f.setRightBrotherIndex(p_node,rb_node);
                else
                    f.setRightBrotherIndex(p_node,0);
                return rbs+1;
            }
        } else {
            h_score = alg_->replace(this->f1_->label(p.i),
                                    getMtrxVal(this->f1_->down(p.i),this->f2_->down(p.k)),
                                    this->f2_->label(p.k),
                                    getMtrxVal(this->f1_->over(p.i,p.j),this->f2_->over(p.k,p.l)));

						//std::cout << "score rep " << h_score << std::endl;
            if (score == h_score) {
								int anchor = 0;
								if (this->anchored_) {
									if (this->f1_->getAnchor(p.i)==this->f2_->getAnchor(p.k)) {
										anchor = this->f1_->getAnchor(p.i);
										//if (anchor != 0) {
										//	std::cout << "found an anchor" << std::endl;
										//	exit(0);
										//}
									}
									else {
										std::cout << "apparently matchend unequal anchors." << std::endl;
										exit(EXIT_FAILURE);
									}
								}

                // it is a replacement
                f.makeRepLabel(p_node,this->f1_->label(p.i),this->f2_->label(p.k),0);
                // down alignment
                noc = backtrace(f, CSFPair(p.i+1,this->f1_->noc(p.i),p.k+1,this->f2_->noc(p.k)), node);
                f.setNumChildren(p_node,noc);
                rb_node = node;
                // right alignment
                rbs = backtrace(f, CSFPair(this->f1_->rb(p.i),p.j-1,this->f2_->rb(p.k),p.l-1), node);
                if (rbs)
                    f.setRightBrotherIndex(p_node,rb_node);
                else
                    f.setRightBrotherIndex(p_node,0);
                return rbs+1;
            }
        }
    }

    // could it be a deletion
    if (p.j>0) {
        h = p.k;               // h is the node where the suffix of the split begins
        for (r=0; r<=p.l; r++) { // for all splits of f2_
            h_score = alg_->del(this->f1_->label(p.i),
                                getMtrxVal(this->f1_->down(p.i),this->f2_->indexpos(p.k,r)),
                                getMtrxVal(this->f1_->over(p.i,p.j),this->f2_->indexpos(h,p.l-r)));

						//std::cout << "score del " << h_score << std::endl;
            if (score == h_score) {
                // it is a deletion
                f.makeDelLabel(p_node,this->f1_->label(p.i));
                // down alignment
                noc = backtrace(f, CSFPair(p.i+1,this->f1_->noc(p.i),p.k,r), node);
                f.setNumChildren(p_node,noc);
                rb_node = node;
                // right alignment
                rbs = backtrace(f, CSFPair(this->f1_->rb(p.i),p.j-1,h,p.l-r), node);
                if (rbs)
                    f.setRightBrotherIndex(p_node,rb_node);
                else
                    f.setRightBrotherIndex(p_node,0);
                return rbs+1;
            }
            //	  if(r<l)       // do not calculate rightbrother of h=0=no-rightbrother
            h = this->f2_->rb(h);
        }
    }

    // could it be an insertion
    if (p.l>0) {
        h = p.i;
        for (r=0; r<=p.j; r++) { // for all splits of f1_
            h_score = alg_->insert(getMtrxVal(this->f1_->indexpos(p.i,r),this->f2_->down(p.k)),
                                   this->f2_->label(p.k),
                                   getMtrxVal(this->f1_->indexpos(h,p.j-r),this->f2_->over(p.k,p.l)));
						//std::cout << "score ins " << h_score << std::endl;
            if (score == h_score) {
                // it is an insertion
                f.makeInsLabel(p_node,this->f2_->label(p.k));
                // down alignment
                noc = backtrace(f, CSFPair(p.i,r,p.k+1,this->f2_->noc(p.k)), node);
                f.setNumChildren(p_node,noc);
                rb_node = node;
                // right alignment
                rbs = backtrace(f, CSFPair(h,p.j-r,this->f2_->rb(p.k),p.l-1), node);
                if (rbs)
                    f.setRightBrotherIndex(p_node,rb_node);
                else
                    f.setRightBrotherIndex(p_node,0);
                return rbs+1;
            }
            //	  if(r<j)
            h = this->f1_->rb(h);
        }
    }

    // you should never be here
    std::cerr << "Strange things happening in backtrace" << std::endl;
		// std::cout << f << std::endl;
    exit(EXIT_FAILURE);
}


// Public functions

template<class R,class L,class AL>
R AlignmentLinear<R,L,AL>::getGlobalOptimum() const {
    return getMtrxVal(this->f1_->indexpos(0,this->f1_->getMaxLength(0)),
											this->f2_->indexpos(0,this->f2_->getMaxLength(0)));
};

template<class R,class L,class AL>
double AlignmentLinear<R,L,AL>::getGlobalOptimumRelative() const {

    double opt = (double) this->getGlobalOptimum();

    double max_x,max_y;
    if (rnaAlg_) {
        max_x = (double) this->f1_->maxScore(*rnaAlg_);
        max_y = (double) this->f2_->maxScore(*rnaAlg_);
    } else {
        max_x = (double) this->f1_->maxScore(*alg_);
        max_y = (double) this->f2_->maxScore(*alg_);
    }

    assert(max_x+max_y>0);
    opt = 2*opt/(max_x+max_y);

    return opt;
};

template<class R,class L,class AL>
R AlignmentLinear<R,L,AL>::getSILOptimum() const {
    R silOptimum;

    if (rnaAlg_)
        silOptimum = rnaAlg_->worst_score();
    else
        silOptimum = alg_->worst_score();

    // find the best match
    unsigned int j = this->f1_->getMaxLength(0);
    unsigned long n = this->f2_->size();
    for (int k=n-1; k>=0; k--)  // for all nodes in f1_
        for (unsigned int l=1; l<=this->f2_->getMaxLength(k); l++) { // for all non empty csfs induced by k
            silOptimum = alg_->choice(silOptimum,getMtrxVal(this->f1_->indexpos(0,j),this->f2_->indexpos(k,l)));
        }

    return silOptimum;
}

template<class R,class L,class AL>
void AlignmentLinear<R,L,AL>::getOptGlobalAlignment(ForestAli<L,AL> &fali) {
    unsigned int node = 0;

    // allocate a forest of the maximal size that a forest alignment can have
    fali.initialize(this->f1_->size()+this->f2_->size(),this->anchored_);
    backtrace(fali,CSFPair(0,this->f1_->getMaxLength(0),0,this->f2_->getMaxLength(0)),node);
    fali.setSize(node);
    fali.calcSumUpCSF();
    fali.calcRMB();
}

template<class R,class L,class AL>
void Alignment<R,L,AL>::resetOptLocalAlignment(int suboptimalsPercent) {
    localAlis_.clear();
    suboptimalsPercent_ = suboptimalsPercent/100.0;
    localSubOptimum_ = localOptimum_;
};

// calculate the score of the next best local alignment that is not "included"
// in a local alignment returned by getOptLocalAlignment before
template<class R,class L,class AL>
bool AlignmentLinear<R,L,AL>::nextLocalSuboptimum() {

    this->localSubOptimum_ = alg_->worst_score();

    // find a matrix element that is optimal
    unsigned long m = this->f1_->size();
    unsigned long n = this->f2_->size();
    for (int i=m-1; i>=0; i--)  // for all nodes in f1_
        for (int k=n-1; k>=0; k--)  // for all nodes in f1_
            for (unsigned int j=1; j<=this->f1_->getMaxLength(i); j++) // for all non empty csfs induced by i
                for (unsigned int l=1; l<=this->f2_->getMaxLength(k); l++) { // for all non empty csfs induced by k
                    bool disjoint=true;

                    // check if i,j,k,l is included
                    typename std::vector<CSFPair>::const_iterator it;
                    for (it=this->localAlis_.begin(); it!=this->localAlis_.end(); it++) {
                        if (!this->f1_->isDisjoint(it->i,it->j,i,j) || !this->f2_->isDisjoint(k,l,it->k,it->l)) {
                            disjoint=false;
                            break;
                        }

                    }

                    if (disjoint) {
                        if (getMtrxVal(this->f1_->indexpos(i,j),this->f2_->indexpos(k,l))>=this->suboptimalsPercent_*this->localOptimum_)
                            this->localSubOptimum_ = alg_->choice(this->localSubOptimum_,getMtrxVal(this->f1_->indexpos(i,j),this->f2_->indexpos(k,l)));
                    }
                }

    return this->localSubOptimum_!=alg_->worst_score();
}

template<class R,class L,class AL>
void AlignmentLinear<R,L,AL>::getOptLocalAlignment(ForestAli<L,AL> &fali,unsigned int &xbasepos, unsigned int &ybasepos) {

    int i = 0, k = 0;
    unsigned int j = 0, l = 0, node = 0;

    // allocate a forest of the maximal size that a forest alignment can have
    fali.initialize(this->f1_->size()+this->f2_->size(),this->anchored_);

    // find a matrix element that is optimal
    unsigned long m = this->f1_->size();
    unsigned long n = this->f2_->size();
    for (i=m-1; i>=0; i--)  // for all nodes in f1_
        for (k=n-1; k>=0; k--)  // for all nodes in f1_
            for (j=1; j<=this->f1_->getMaxLength(i); j++) // for all non empty csfs induced by i
                for (l=1; l<=this->f2_->getMaxLength(k); l++) { // for all non empty csfs induced by k
                    bool disjoint = true;

                    // check if i,j,k,l is included
                    typename std::vector<CSFPair>::const_iterator it;
                    for (it=this->localAlis_.begin(); it!=this->localAlis_.end(); it++) {
                        if (!this->f1_->isDisjoint(it->i,it->j,i,j) || !this->f2_->isDisjoint(k,l,it->k,it->l)) {
                            disjoint = false;
                            break;
                        }
                    }
                    if (disjoint && getMtrxVal(this->f1_->indexpos(i,j),this->f2_->indexpos(k,l)) == this->localSubOptimum_)
                        goto found;
                }

found:
    backtrace(fali,CSFPair(i,j,k,l),node);
    fali.setSize(node);
    fali.calcSumUpCSF();
    fali.calcRMB();

    this->localAlis_.push_back(CSFPair(i,j,k,l));
    xbasepos=this->f1_->countLeaves(i);
    ybasepos=this->f2_->countLeaves(k);
}

template<class R,class L,class AL>
void AlignmentLinear<R,L,AL>::getOptSILAlignment(ForestAli<L,AL> &fali,unsigned int &ybasepos) {
    unsigned int node=0;

    // allocate a forest of the maximal size that a forest alignment can have
    fali.initialize(this->f1_->size()+this->f2_->size(),this->anchored_);

    R silOptimum = getSILOptimum();
    unsigned int j = this->f1_->getMaxLength(0);
    int k = 0;
    unsigned int l = 0;

    // find a matrix element that is optimal
    unsigned long n = this->f2_->size();
    for (k=n-1; k>=0; k--)  // for all nodes in f1_
        for (l=1; l<=this->f2_->getMaxLength(k); l++) { // for all non empty csfs induced by k
            if (silOptimum == getMtrxVal(this->f1_->indexpos(0,j),this->f2_->indexpos(k,l)))
                goto found;
        }

found:
    backtrace(fali,CSFPair(0,j,k,l),node);
    fali.setSize(node);
    fali.calcSumUpCSF();
    fali.calcRMB();

    ybasepos = this->f2_->countLeaves(k);
}


template<class R, class L, class AL>
AlignmentAffine<R,L,AL>::AlignmentAffine(const Forest<L> *f1, const Forest<L> *f2, const AlgebraAffine<R,L> &alg, 
		const bool topdown, const bool anchored, const bool local, const bool printBacktrace, bool speedup)
        : Alignment<R,L,AL>(f1, f2, topdown, anchored, printBacktrace) {

    // alloc space for the score matrix, backtrace structure,
    // and , if wanted, for the calculation-order-matrix
		mtrx_ = new TAD_DP_TableAffine<R>(this->f1_->getNumCSFs(),this->f2_->getNumCSFs(),alg.worst_score());
    // initialize variables
    alg_ = &alg;
    rnaAlg_ = NULL;
    this->localOptimum_ = alg.worst_score();
    // align forests f1_ and f2_
		bool RNA = false;
		if (topdown)
			this->fillMatricesTopDown(local, RNA, speedup);
		else
			this->fillMatricesBottomUp(local, RNA, speedup);
}


// constructor for RNA alignments
template<class R,class L,class AL>
AlignmentAffine<R,L,AL>::AlignmentAffine(const Forest<L> *f1, const Forest<L> *f2, const RNA_AlgebraAffine<R,L> &rnaAlg, 
		const bool topdown, const bool anchored, const bool local, const bool printBacktrace, bool speedup)
        : Alignment<R,L,AL>(f1, f2, topdown, anchored, printBacktrace) {

    // alloc space for the score matrix, backtrace structure and,
    // if wanted, for the calculation-order-matrix
		mtrx_ = new TAD_DP_TableAffine<R>(this->f1_->getNumCSFs(),this->f2_->getNumCSFs(),rnaAlg.worst_score());
    // initialize variables
    rnaAlg_ = &rnaAlg;
    alg_ = (const AlgebraAffine<R,L>*)&rnaAlg;
    this->localOptimum_ = rnaAlg.worst_score();
    // align forests f1 and f2
		bool RNA = true;
		if (topdown)
			this->fillMatricesTopDown(local, RNA, speedup);
		else
			this->fillMatricesBottomUp(local, RNA, speedup);
}

template<class R, class L, class AL>
void AlignmentAffine<R, L, AL>::makeFirstCell() {
  // the easiest case ..
  for (int t = S; t <= VH_; t++) // for all tables
    setMtrxVal(t, 0, 0, alg_->empty());
	//std::cout << "affine!!" << std::endl;
}

template<class R, class L, class AL>
void AlignmentAffine<R, L, AL>::makeFirstCol() {

  // align f1 to the empty forest (fill first column of array)
  // for all nodes in f1
  unsigned long m = this->f1_->size();
  for (long i = m - 1; i >= 0; i--) {
    // for all non empty csfs induced by i
    for (unsigned long j = 1; j <= this->f1_->getMaxLength(i); j++) {

			// temporary info for next tables and opening or extend
			unsigned int t_c = 0;
			unsigned int t_rb = 0;
			bool open = false;
			R score = alg_->worst_score(); 

      for (int t = S; t <= VH_; t++) {// for all tables
				getConnectedTablesDel(t, t_c, t_rb, open);
				R dwn = getMtrxVal(t_c, this->f1_->down(i), 0);
				R ovr = getMtrxVal(t_rb, this->f1_->over(i,j), 0);
				L label_i = this->f1_->label(i);
				if (open)
					score = alg_->delO(label_i, dwn, ovr);
				else
					score = alg_->del(label_i, dwn, ovr);
        setMtrxVal(t,this->f1_->indexpos(i, j), 0, score);
			}

    }
  }

}

template<class R, class L, class AL>
void AlignmentAffine<R, L, AL>::makeFirstRow() {

  // align f2 to the empty forest (fill first row of array)
  // for all nodes in f1
  unsigned long n = this->f2_->size();
  for (long k = n - 1; k >= 0; k--) {
    // for all non empty csfs induced by k
    for (unsigned long l = 1; l <= this->f2_->getMaxLength(k); l++) {

			// temporary info for next tables and opening or extend
			unsigned int t_c = 0;
			unsigned int t_rb = 0;
			bool open = false;
			R score = alg_->worst_score(); 

      for (int t = S; t <= VH_; t++) {// for all tables
				getConnectedTablesIns(t, t_c, t_rb, open);
				R dwn = getMtrxVal(t_c, 0, this->f2_->down(k));
				R ovr = getMtrxVal(t_rb, 0, this->f2_->over(k, l));
				L label_k = this->f2_->label(k);
				if (open)
					score = alg_->insertO(dwn, label_k, ovr);
				else
					score = alg_->insert(dwn, label_k, ovr);
        setMtrxVal(t, 0, this->f2_->indexpos(k, l), score);
			}

    }
  }

}


template<class R, class L, class AL>
void AlignmentAffine<R, L, AL>::recursiveFillAnchored(CSFPair p, bool RNA, bool speedup) {
		//std::cout << "rec fill anch not imp yet" << std::endl;
    std::vector<R> score(7, this->alg_->worst_score());
    std::vector<R> h_score(7, this->alg_->worst_score());
    int a = this->f1_->getAnchor(p.i);
    int b = this->f2_->getAnchor(p.k);
    //this->calls++;

    // easiest: already done;
    if (computed(this->f1_->indexpos(p.i,p.j) , this->f2_->indexpos(p.k, p.l)))
        return;

    // empty alignment
    // no recursion for empty ali
    if (p.j == 0 && p.l == 0) {
        this->makeFirstCell();
        return;
    }

    // replacement
    // a and b both anchors or both nonanchors
    if (a==b && p.j > 0 && p.l > 0) {
        if (a>0 && b>0) {
            //std::cout << "Matched anchor " << a << "." << std::endl;
        }
    		if (RNA && this->f1_->down(p.i) && this->f2_->down(p.k)) {
						if (! computed(this->f1_->mdown(p.i),this->f2_->mdown(p.k)))
							recursiveFillAnchored(CSFPair(p.i + 1 + 1, this->f1_->noc(p.i) - 2, p.k + 1 + 1, this->f2_->noc(p.k) - 2), RNA, speedup);
        }
				else {
						if (! computed(this->f1_->down(p.i), this->f2_->down(p.k)))
							recursiveFillAnchored(CSFPair(p.i + 1, this->f1_->noc(p.i), p.k + 1, this->f2_->noc(p.k)), RNA, speedup);
				}
				if (! computed(this->f1_->over(p.i, p.j) , this->f2_->over(p.k, p.l)))
					recursiveFillAnchored(CSFPair(this->f1_->rb(p.i), p.j - 1, this->f2_->rb(p.k), p.l - 1), RNA, speedup);
				// TODO here we now immediately compute the replacement..

			  for (int t = S; t <= VH_; t++) { // for all tables
		
		      unsigned int t_c = 0;
		      unsigned int t_rb = 0;
		      getConnectedTablesRep(t, t_c, t_rb);
		
					R ovr = getMtrxVal(t_rb, this->f1_->over(p.i,p.j),this->f2_->over(p.k,p.l));
					// basepair replacement
				  // must be two P nodes !!
					if (rnaAlg_ && this->f1_->down(p.i) && this->f2_->down(p.k)) {
						R dwn = getMtrxVal(t_c, this->f1_->mdown(p.i), this->f2_->mdown(p.k));
			      score[t] = rnaAlg_->replacepair(this->f1_->label(p.i+1), this->f2_->label(p.k+1), dwn,
			                   this->f1_->label(this->f1_->getRightmostBrotherIndex(p.i+1)),
			                   this->f2_->label(this->f2_->getRightmostBrotherIndex(p.k+1)), ovr);
					}// base replacement
					else {
						R dwn = getMtrxVal(t_c, this->f1_->down(p.i),this->f2_->down(p.k));
			      score[t] = alg_->replace(this->f1_->label(p.i), dwn, this->f2_->label(p.k), ovr);
			    }
				}

				// replacement ende
				// two anchors were matched - done
        if (a>0 && b>0) {
          // set value
          for (int t = S; t <= VH_; t++) // for all tables
            setMtrxVal(t,this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l), score[t]);
          setComputed(this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l));
          return;
        }
		}
		
    // deletion
    if (a==0 && p.j > 0) {
	      for (int t = S; t <= VH_; t++) { // for all tables
	
	        bool open = false;
	        unsigned int t_c = 0;
	        unsigned int t_rb = 0;
	        getConnectedTablesDel(t, t_c, t_rb, open);
	
					L lbl = this->f1_->label(p.i);
	
			    if (this->f1_->noc(p.i)==0 && speedup) { // no child
							// right alignment
							if (! computed(this->f1_->over(p.i, p.j) , this->f2_->indexpos(p.k, p.l)))
								recursiveFillAnchored(CSFPair(this->f1_->rb(p.i), p.j - 1, p.k, p.l), RNA, speedup);
							R ovr = getMtrxVal(t_rb, this->f1_->over(p.i,p.j),this->f2_->indexpos(p.k,p.l));
							if (open)
								h_score[t] = alg_->delO(lbl, 0, ovr);
							else
								h_score[t] = alg_->del(lbl, 0, ovr);
			        score[t] = alg_->choice(score[t],h_score[t]);
			    } 
					else  if (this->f1_->rb(p.i)==0 && speedup) { // no right brother
							// down alignment
							if (! computed(this->f1_->down(p.i) , this->f2_->indexpos(p.k, p.l)))
								recursiveFillAnchored(CSFPair(p.i + 1, this->f1_->noc(p.i), p.k, p.l), RNA, speedup);
							R dwn = getMtrxVal(t_c, this->f1_->down(p.i),this->f2_->indexpos(p.k,p.l));
							if (open)
								h_score[t] = alg_->delO(lbl, dwn, 0);
							else
								h_score[t] = alg_->del(lbl, dwn, 0);
							score[t] = alg_->choice(score[t],h_score[t]);
			    } 
					else {
						// h = node where suffix of split begins
						unsigned int h = p.k;
						// for all splits of f2
				  	for (unsigned int r = 0; r <= p.l; r++) {
								// down alignment
								if (! computed(this->f1_->down(p.i) , this->f2_->indexpos(p.k, r)))
									recursiveFillAnchored(CSFPair(p.i + 1, this->f1_->noc(p.i), p.k, r), RNA, speedup);
								// right alignment
								if (! computed(this->f1_->over(p.i, p.j) , this->f2_->indexpos(h, p.l - r)))
									recursiveFillAnchored(CSFPair(this->f1_->rb(p.i), p.j - 1, h, p.l - r), RNA, speedup);
								R dwn = getMtrxVal(t_c, this->f1_->down(p.i), this->f2_->indexpos(p.k, r));
								R ovr = getMtrxVal(t_rb, this->f1_->over(p.i, p.j), this->f2_->indexpos(h, p.l - r));
						    if (open) 
						        h_score[t] = alg_->delO(lbl, dwn, ovr);
						    else 
						        h_score[t] = alg_->del(lbl, dwn, ovr);
						      
				  	    score[t] = alg_->choice(score[t], h_score[t]);
								h = this->f2_->rb(h);
						}
					}
				}
    }

    // insertion
    if (b==0 && p.l > 0) {
	      for (int t = S; t <= VH_; t++) { // for all tables
	
	        bool open = false;
	        unsigned int t_c = 0;
	        unsigned int t_rb = 0;
	        getConnectedTablesIns(t, t_c, t_rb, open);
	
					L lbl = this->f2_->label(p.k);
			    if (this->f2_->noc(p.k)==0 && speedup) { // no child
							// right alignment
							if (! computed(this->f1_->indexpos(p.i, p.j) , this->f2_->over(p.k, p.l)))
								recursiveFillAnchored(CSFPair(p.i, p.j, this->f2_->rb(p.k), p.l - 1), RNA, speedup);
							R ovr = getMtrxVal(t_rb, this->f1_->indexpos(p.i,p.j),this->f2_->over(p.k,p.l));
							if (open)
								h_score[t] = alg_->insertO(0, lbl, ovr);
							else 
								h_score[t] = alg_->insert(0, lbl, ovr);
			        score[t] = alg_->choice(score[t],h_score[t]);
			    } 
					else if (this->f2_->rb(p.k)==0 && speedup) { // no right brother
							// down alignment
							if (! computed(this->f1_->indexpos(p.i, p.j) , this->f2_->down(p.k)))
								recursiveFillAnchored(CSFPair(p.i, p.j, p.k + 1, this->f2_->noc(p.k)), RNA, speedup);
							R dwn = getMtrxVal(t_c, this->f1_->indexpos(p.i,p.j),this->f2_->down(p.k));
						  if (open) 
								h_score[t] = alg_->insertO(dwn, lbl, 0);
							else
								h_score[t] = alg_->insert(dwn, lbl, 0);
			        score[t] = alg_->choice(score[t],h_score[t]);
			    } 
					else {
							// h = node where suffix of split begins
						  unsigned int h = p.i;
						  // for all splits of f1
						  for (unsigned int r = 0; r <= p.j; r++) {
									// down alignment
									if (! computed(this->f1_->indexpos(p.i, r) , this->f2_->down(p.k)))
										recursiveFillAnchored(CSFPair(p.i, r, p.k + 1, this->f2_->noc(p.k)), RNA, speedup);
									// right alignment
									if (! computed(this->f1_->indexpos(h, p.j - r) , this->f2_->over(p.k, p.l)))
										recursiveFillAnchored(CSFPair(h, p.j - r, this->f2_->rb(p.k), p.l - 1), RNA, speedup);
									R dwn = getMtrxVal(t_c,this->f1_->indexpos(p.i, r), this->f2_->down(p.k));
									R ovr = getMtrxVal(t_rb,this->f1_->indexpos(h, p.j - r), this->f2_->over(p.k, p.l));
						      if (open) 
						      	  h_score[t] = alg_->insertO(dwn, lbl, ovr);
						      else 
						      	  h_score[t] =  alg_->insert(dwn, lbl, ovr);
									//std::cout << "dwn " << dwn << " ovr " << ovr << " open " << open << std::endl;
									//std::cout << "score t " << score[t] << " h_score t " << h_score[t] << std::endl;
						      score[t] = alg_->choice(score[t], h_score[t]);
									h = this->f1_->rb(h);
						  }
			    }
	      }
    }
    // set value
    for (int t = S; t <= VH_; t++) {// for all tables
			if (score[t] > 147483580) {
				std::cout << "score overflow problem - score " << std::endl;
				std::cout << score[t] << " seen in matrix " << t << std::endl;
				exit(0);
			}
      setMtrxVal(t,this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l), score[t]);
		}

    setComputed(this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l));

}


// TODO speedup als parameter
template<class R, class L, class AL>
void AlignmentAffine<R, L, AL>::compareCSFPairTD(CSFPair p, bool speedup) {
  std::vector<R> score(7, alg_->worst_score());
  std::vector<R> h_score(7, alg_->worst_score());

	//std::cout << "topdown" << std::endl;
  // empty alignment
  if (p.j == 0 && p.l == 0) {
    makeFirstCell();
		setComputed(this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l));
    return;
  }
  // replacement
  if (p.j > 0 && p.l > 0) {
	  for (int t = S; t <= VH_; t++) { // for all tables
	
	      unsigned int t_c = 0;
	      unsigned int t_rb = 0;
	      getConnectedTablesRep(t, t_c, t_rb);
	
				R ovr = getMtrxVal(t_rb, this->f1_->over(p.i,p.j),this->f2_->over(p.k,p.l));
				// basepair replacement
			  // must be two P nodes !!
				if (rnaAlg_ && this->f1_->down(p.i) && this->f2_->down(p.k)) {
					R dwn = getMtrxVal(t_c, this->f1_->mdown(p.i), this->f2_->mdown(p.k));
		      score[t] = rnaAlg_->replacepair(this->f1_->label(p.i+1), this->f2_->label(p.k+1), dwn,
		                   this->f1_->label(this->f1_->getRightmostBrotherIndex(p.i+1)),
		                   this->f2_->label(this->f2_->getRightmostBrotherIndex(p.k+1)), ovr);
				}// base replacement
				else {
					R dwn = getMtrxVal(t_c, this->f1_->down(p.i),this->f2_->down(p.k));
		      score[t] = alg_->replace(this->f1_->label(p.i), dwn, this->f2_->label(p.k), ovr);
		    }
	  }
  }
  // delete
  if (p.j > 0) {
    // for all splits of f2
    for (unsigned long r = 0; r <= p.l; r++) {
      for (int t = S; t <= VH_; t++) { // for all tables

        bool open = false;
        unsigned int t_c = 0;
        unsigned int t_rb = 0;
        getConnectedTablesDel(t, t_c, t_rb, open);

				L lbl = this->f1_->label(p.i);

		    if (this->f1_->noc(p.i)==0 && speedup) { // no child
						R ovr = getMtrxVal(t_rb, this->f1_->over(p.i,p.j),this->f2_->indexpos(p.k,p.l));
						if (open)
							h_score[t] = alg_->delO(lbl, 0, ovr);
						else
							h_score[t] = alg_->del(lbl, 0, ovr);
		        score[t] = alg_->choice(score[t],h_score[t]);
		    } 
				else  if (this->f1_->rb(p.i)==0 && speedup) { // no right brother
						R dwn = getMtrxVal(t_c, this->f1_->down(p.i),this->f2_->indexpos(p.k,p.l));
						if (open)
							h_score[t] = alg_->delO(lbl, dwn, 0);
						else
							h_score[t] = alg_->del(lbl, dwn, 0);
						score[t] = alg_->choice(score[t],h_score[t]);
		    } 
				else {
					// h = node where suffix of split begins
					unsigned int h = p.k;
					// for all splits of f2
			  	for (unsigned int r = 0; r <= p.l; r++) {
							R dwn = getMtrxVal(t_c, this->f1_->down(p.i), this->f2_->indexpos(p.k, r));
							R ovr = getMtrxVal(t_rb, this->f1_->over(p.i, p.j), this->f2_->indexpos(h, p.l - r));
					    if (open) 
					        h_score[t] = alg_->delO(lbl, dwn, ovr);
					    else 
					        h_score[t] = alg_->del(lbl, dwn, ovr);
					      
			  	    score[t] = alg_->choice(score[t], h_score[t]);
							h = this->f2_->rb(h);
					}
				}
			}
		}
  }
  // insert
  if (p.l > 0) {
    // for all splits of f1
    for (unsigned long r = 0; r <= p.j; r++) {
      for (int t = S; t <= VH_; t++) { // for all tables

        bool open = false;
        unsigned int t_c = 0;
        unsigned int t_rb = 0;
        getConnectedTablesIns(t, t_c, t_rb, open);

				L lbl = this->f2_->label(p.k);
		    if (this->f2_->noc(p.k)==0 && speedup) { // no child
						R ovr = getMtrxVal(t_rb, this->f1_->indexpos(p.i,p.j),this->f2_->over(p.k,p.l));
						if (open)
							h_score[t] = alg_->insertO(0, lbl, ovr);
						else 
							h_score[t] = alg_->insert(0, lbl, ovr);
		        score[t] = alg_->choice(score[t],h_score[t]);
		    } 
				else if (this->f2_->rb(p.k)==0 && speedup) { // no right brother
						R dwn = getMtrxVal(t_c, this->f1_->indexpos(p.i,p.j),this->f2_->down(p.k));
					  if (open) 
							h_score[t] = alg_->insertO(dwn, lbl, 0);
						else
							h_score[t] = alg_->insert(dwn, lbl, 0);
		        score[t] = alg_->choice(score[t],h_score[t]);
		    } 
				else {
						// h = node where suffix of split begins
					  unsigned int h = p.i;
					  // for all splits of f1
					  for (unsigned int r = 0; r <= p.j; r++) {
								R dwn = getMtrxVal(t_c,this->f1_->indexpos(p.i, r), this->f2_->down(p.k));
								R ovr = getMtrxVal(t_rb,this->f1_->indexpos(h, p.j - r), this->f2_->over(p.k, p.l));
					      if (open) 
					      	  h_score[t] = alg_->insertO(dwn, lbl, ovr);
					      else 
					      	  h_score[t] =  alg_->insert(dwn, lbl, ovr);
					      score[t] = alg_->choice(score[t], h_score[t]);
								h = this->f1_->rb(h);
					  }
		    }
      }
    }
  }

  // set value
  for (int t = S; t <= VH_; t++) // for all tables
    setMtrxVal(t,this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l), score[t]);

  setComputed(this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l));
}

template<class R, class L, class AL>
void AlignmentAffine<R, L, AL>::compareCSFPair(CSFPair p, bool RNA, bool speedup) {

  for (int t = S; t <= VH_; t++) { // for all tables

		// temporary info for next tables and opening or extend
		unsigned int t_c = 0;
		unsigned int t_rb = 0;
		bool open = false;

		R score = alg_->worst_score();
		R h_score = alg_->worst_score();

		// replacement
		getConnectedTablesRep(t, t_c, t_rb);
		R ovr = getMtrxVal(t_rb, this->f1_->over(p.i,p.j),this->f2_->over(p.k,p.l));
		// basepair replacement
	  // must be two P nodes !!
		if (rnaAlg_ && this->f1_->down(p.i) && this->f2_->down(p.k)) {
			R dwn = getMtrxVal(t_c, this->f1_->mdown(p.i), this->f2_->mdown(p.k));
      score = rnaAlg_->replacepair(this->f1_->label(p.i+1), this->f2_->label(p.k+1), dwn,
                   this->f1_->label(this->f1_->getRightmostBrotherIndex(p.i+1)),
                   this->f2_->label(this->f2_->getRightmostBrotherIndex(p.k+1)), ovr);
		}// base replacement
		else {
			R dwn = getMtrxVal(t_c, this->f1_->down(p.i),this->f2_->down(p.k));
      score = alg_->replace(this->f1_->label(p.i), dwn, this->f2_->label(p.k), ovr);
    }

	  // h = node where suffix of split begins
	  getConnectedTablesDel(t, t_c, t_rb, open);
		L lbl = this->f1_->label(p.i);

    if (this->f1_->noc(p.i)==0 && speedup) { // no child
				R ovr = getMtrxVal(t_rb, this->f1_->over(p.i,p.j),this->f2_->indexpos(p.k,p.l));
				if (open)
					h_score = alg_->delO(lbl, 0, ovr);
				else
					h_score = alg_->del(lbl, 0, ovr);
        score = alg_->choice(score,h_score);
    } 
		else  if (this->f1_->rb(p.i)==0 && speedup) { // no right brother
				R dwn = getMtrxVal(t_c, this->f1_->down(p.i),this->f2_->indexpos(p.k,p.l));
				if (open)
					h_score = alg_->delO(lbl, dwn, 0);
				else
					h_score = alg_->del(lbl, dwn, 0);
				score = alg_->choice(score,h_score);
    } 
		else {
			unsigned int h = p.k;
			// for all splits of f2
	  	for (unsigned int r = 0; r <= p.l; r++) {
					R dwn = getMtrxVal(t_c, this->f1_->down(p.i), this->f2_->indexpos(p.k, r));
					R ovr = getMtrxVal(t_rb, this->f1_->over(p.i, p.j), this->f2_->indexpos(h, p.l - r));
			    if (open) 
			        h_score = alg_->delO(lbl, dwn, ovr);
			    else 
			        h_score = alg_->del(lbl, dwn, ovr);
			      
	  	    score = alg_->choice(score, h_score);
					h = this->f2_->rb(h);
	  	}
    } 

	  // insert
	  getConnectedTablesIns(t, t_c, t_rb, open);
		lbl = this->f2_->label(p.k);
    if (this->f2_->noc(p.k)==0 && speedup) { // no child
				R ovr = getMtrxVal(t_rb, this->f1_->indexpos(p.i,p.j),this->f2_->over(p.k,p.l));
				if (open)
					h_score = alg_->insertO(0, lbl, ovr);
				else 
					h_score = alg_->insert(0, lbl, ovr);
        score = alg_->choice(score,h_score);
    } 
		else if (this->f2_->rb(p.k)==0 && speedup) { // no right brother
				R dwn = getMtrxVal(t_c, this->f1_->indexpos(p.i,p.j),this->f2_->down(p.k));
			  if (open) 
					h_score = alg_->insertO(dwn, lbl, 0);
				else
					h_score = alg_->insert(dwn, lbl, 0);
        score = alg_->choice(score,h_score);
    } 
		else {
			  unsigned int h = p.i;
			  // for all splits of f1
			  for (unsigned int r = 0; r <= p.j; r++) {
						R dwn = getMtrxVal(t_c,this->f1_->indexpos(p.i, r), this->f2_->down(p.k));
						R ovr = getMtrxVal(t_rb,this->f1_->indexpos(h, p.j - r), this->f2_->over(p.k, p.l));
			      if (open) 
			      	  h_score = alg_->insertO(dwn, lbl, ovr);
			      else 
			      	  h_score =  alg_->insert(dwn, lbl, ovr);
			      score = alg_->choice(score, h_score);
						h = this->f1_->rb(h);
			  }
    }

	  // set value
    setMtrxVal(t,this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l), score);
	}
}


// to which table do we go next in backtracing?
// do we open or extend the gaps?
template<class R, class L, class AL>
void AlignmentAffine<R, L, AL>::getConnectedTablesRep(const unsigned int t,
    unsigned int &t_c, unsigned int &t_rb) {

  t_c = S;
  switch (t) {
    case S: t_rb = S; break;
    case V: t_rb = V; break;
    case H: t_rb = H; break;
    case V_: t_rb = S; break;
    case H_: t_rb = S; break;
    case V_H: t_rb = H; break;
    case VH_: t_rb = V; break;
    default:
      std::cerr << t << " is no valid table for backtracing." << std::endl;
      exit(EXIT_FAILURE);
  }
}

// to which table do we go next in backtracing?
// do we open or extend the gaps?
template<class R, class L, class AL>
void AlignmentAffine<R, L, AL>::getConnectedTablesDel(const unsigned int t,
    unsigned int &t_c, unsigned int &t_rb, bool &open) {
  t_c= H;
  open= false;
  switch (t) {
    case S: t_rb = H_; open = true; break;
    case V: t_rb = VH_; open = true; break;
    case H: t_rb = H; break;
    case V_: t_rb = H_; open = true; break;
    case H_: t_rb = H_; break;
    case V_H: t_rb = H; break;
    case VH_: t_rb = VH_; break;
    default:
      std::cerr << t << " is no valid table for backtracing." << std::endl;
      exit(EXIT_FAILURE);
  }
}

// to which table do we go next in backtracing?
// do we open or extend the gaps?
template<class R, class L, class AL>
void AlignmentAffine<R, L, AL>::getConnectedTablesIns(const unsigned int t,
    unsigned int &t_c, unsigned int &t_rb, bool &open) {
  t_c= V;
  open= false;
  switch (t) {
    case S: t_rb = V_; open = true; break;
    case V: t_rb = V; break;
    case H: t_rb = V_H; open = true; break;
    case V_: t_rb = V_; break;
    case H_: t_rb = V_; open = true; break;
    case V_H: t_rb = V_H; break;
    case VH_: t_rb = V; break;
    default:
      std::cerr << t << " is no valid table for backtracing." << std::endl;
      exit(EXIT_FAILURE);
  }
}

template<class R, class L, class AL>
R AlignmentAffine<R, L, AL>::getGlobalOptimum() const {
	  return getMtrxVal(S, this->f1_->indexpos(0, this->f1_->getMaxLength(0)),
				                 this->f2_->indexpos(0, this->f2_->getMaxLength(0)));
}

template<class R, class L, class AL>
double AlignmentAffine<R, L, AL>::getGlobalOptimumRelative() const {

  double opt = (double) getGlobalOptimum();

  double max_x,max_y;
  if (rnaAlg_) {
      max_x = (double) this->f1_->maxScore(*rnaAlg_);
      max_y = (double) this->f2_->maxScore(*rnaAlg_);
  } else {
      max_x = (double) this->f1_->maxScore(*alg_);
      max_y = (double) this->f2_->maxScore(*alg_);
  }

  assert(max_x + max_y > 0);
  opt = 2 * opt / (max_x + max_y);

  return opt;
}


template<class R,class L,class AL>
R AlignmentAffine<R,L,AL>::getSILOptimum() const {
    R silOptimum;

    if (rnaAlg_)
        silOptimum = rnaAlg_->worst_score();
    else
        silOptimum = alg_->worst_score();

		// TODO in which table(s) can we start?
    // find the best match
    unsigned int j = this->f1_->getMaxLength(0);
    unsigned long n = this->f2_->size();
    for (int k=n-1; k>=0; k--)  // for all nodes in f1_
        for (unsigned int l=1; l<=this->f2_->getMaxLength(k); l++) { // for all non empty csfs induced by k
            silOptimum = alg_->choice(silOptimum,getMtrxVal(S,this->f1_->indexpos(0,j),this->f2_->indexpos(k,l)));
        }

    return silOptimum;
}


// print the alignment in an alignment forest
template<class R, class L, class AL>
void AlignmentAffine<R, L, AL>::getOptGlobalAlignment(ForestAli<L, AL> &fali) {
  unsigned int node = 0;

  // allocate a forest of the maximal size that a forest alignment can have
  fali.initialize(this->f1_->size() + this->f2_->size(),this->anchored_);
	
  backtrace(fali, CSFPair(0, this->f1_->getMaxLength(0), 0, this->f2_->getMaxLength(0)), node, S);
  fali.setSize(node);
  fali.calcSumUpCSF();
  fali.calcRMB();

}


template<class R, class L, class AL>
void AlignmentAffine<R, L, AL>::getOptLocalAlignment(ForestAli<L, AL> &fali,
    unsigned int &start1, unsigned int &start2) {

    int i = 0, k = 0;
    unsigned int j = 0, l = 0, node = 0;

    // allocate a forest of the maximal size that a forest alignment can have
    fali.initialize(this->f1_->size()+this->f2_->size(),this->anchored_);

    // find a matrix element that is optimal
    unsigned long m = this->f1_->size();
    unsigned long n = this->f2_->size();
    for (i=m-1; i>=0; i--)  // for all nodes in f1_
        for (k=n-1; k>=0; k--)  // for all nodes in f1_
            for (j=1; j<=this->f1_->getMaxLength(i); j++) // for all non empty csfs induced by i
                for (l=1; l<=this->f2_->getMaxLength(k); l++) { // for all non empty csfs induced by k
                    bool disjoint = true;

                    // check if i,j,k,l is included
                    typename std::vector<CSFPair>::const_iterator it;
                    for (it=this->localAlis_.begin(); it!=this->localAlis_.end(); it++) {
                        if (!this->f1_->isDisjoint(it->i,it->j,i,j) || !this->f2_->isDisjoint(k,l,it->k,it->l)) {
                            disjoint = false;
                            break;
                        }
                    }
                    if (disjoint && getMtrxVal(this->localOptimumTable_,this->f1_->indexpos(i,j),this->f2_->indexpos(k,l)) == this->localSubOptimum_)
                        goto found;
                }
found:
    backtrace(fali,CSFPair(i,j,k,l),node, this->localOptimumTable_);
    fali.setSize(node);
    fali.calcSumUpCSF();
    fali.calcRMB();

    this->localAlis_.push_back(CSFPair(i,j,k,l));
    start1 = this->f1_->countLeaves(i);
    start2 = this->f2_->countLeaves(k);
}

template<class R, class L, class AL>
void AlignmentAffine<R, L, AL>::foundOptLocalAlignment(ForestAli<L, AL> &fali, CSFPair p,
    unsigned int &s1, unsigned int &s2,
    unsigned int &e1, unsigned int &e2) {

  unsigned int node = 0;
  //if (options_.printBacktrack) {
  //  std::cout << "backtracing of local affine alignment starting at pos " << this->f1_->indexpos(i, p.j) << ", " << this->f2_->indexpos(p.k, p.l) << std::endl;
  //}

  unsigned int start1 = p.i;
  unsigned int start2 = p.k;
  s1 = this->f1_->countLeaves(start1);
  s2 = this->f2_->countLeaves(start2);

  backtrace(fali, p, node, this->localOptimumTable_);
  fali.setSize(node);
  fali.calcSumUpCSF();
  fali.calcRMB();
  this->localAlis_.push_back(p);

	// TODO from where do we have this info?
  //e1 = this->f1_->countLeaves(end1);
  //e2 = this->f2_->countLeaves(end2);
}


// calculate the score of the next best local alignment that is not "included"
// in a local alignment returned by getOptLocalAlignment before
template<class R,class L,class AL>
bool AlignmentAffine<R,L,AL>::nextLocalSuboptimum() {

    this->localSubOptimum_ = alg_->worst_score();

    // find a matrix element that is optimal
    unsigned long m = this->f1_->size();
    unsigned long n = this->f2_->size();
    for (int i=m-1; i>=0; i--)  // for all nodes in f1_
        for (int k=n-1; k>=0; k--)  // for all nodes in f1_
            for (unsigned int j=1; j<=this->f1_->getMaxLength(i); j++) // for all non empty csfs induced by i
                for (unsigned int l=1; l<=this->f2_->getMaxLength(k); l++) { // for all non empty csfs induced by k
                    bool disjoint=true;

                    // check if i,j,k,l is included
                    typename std::vector<CSFPair>::const_iterator it;
                    for (it=this->localAlis_.begin(); it!=this->localAlis_.end(); it++) {
                        if (!this->f1_->isDisjoint(it->i,it->j,i,j) || !this->f2_->isDisjoint(k,l,it->k,it->l)) {
                            disjoint=false;
                            break;
                        }

                    }

										// TODO we only look at matrix S, so the alignment can only end in "normal" state (no gap state)
                    if (disjoint) {
                        if (getMtrxVal(S,this->f1_->indexpos(i,j),this->f2_->indexpos(k,l))>=this->suboptimalsPercent_*this->localOptimum_)
                            this->localSubOptimum_ = alg_->choice(this->localSubOptimum_,getMtrxVal(S,this->f1_->indexpos(i,j),this->f2_->indexpos(k,l)));
                    }
                }

    return this->localSubOptimum_!=alg_->worst_score();
}


template<class R,class L,class AL>
void AlignmentAffine<R,L,AL>::getOptSILAlignment(ForestAli<L,AL> &fali,unsigned int &ybasepos) {
    unsigned int node=0;

    // allocate a forest of the maximal size that a forest alignment can have
    fali.initialize(this->f1_->size()+this->f2_->size(), this->anchored_);

    R silOptimum = getSILOptimum();
    unsigned int j = this->f1_->getMaxLength(0);
    int k = 0;
    unsigned int l = 0;

    // find a matrix element that is optimal
    unsigned long n = this->f2_->size();
    for (k=n-1; k>=0; k--)  // for all nodes in f1_
        for (l=1; l<=this->f2_->getMaxLength(k); l++) { // for all non empty csfs induced by k
						// TODO we only look at matrix S, so the alignment can only end in "normal" state (no gap state)
            if (silOptimum == getMtrxVal(S,this->f1_->indexpos(0,j),this->f2_->indexpos(k,l)))
                goto found;
        }

found:
    backtrace(fali,CSFPair(0,j,k,l),node,S);
    fali.setSize(node);
    fali.calcSumUpCSF();
    fali.calcRMB();

    ybasepos = this->f2_->countLeaves(k);
}

// NORMAL BACKTRACE, not for RNA

template<class R, class L, class AL>
unsigned int AlignmentAffine<R,L,AL>::backtrace(ForestAli<L,AL> &f, CSFPair p, unsigned int &node, int t) {

		if (this->printBacktrace_)
			std::cout << "backtracing from i " << p.i << " j " << p.j << " k " << p.k << " l " << p.l   << std::endl;

//  if (p.i > end1)
//    end1 = p.i; // keep track of nodes in alignment - we need the last nodes in the end
//  if (p.k > end2)
//    end2 = p.k;

  unsigned int p_node = 0;//, rb_node = 0;
  unsigned int r = 0, h = 0;//noc = 0, rbs = 0;

  // empty AlignmentAffine
  if (p.j == 0 && p.l == 0) {
    return 0;
  }

  R h_score = alg_->worst_score();
  R score = getMtrxVal(t, this->f1_->indexpos(p.i, p.j), this->f2_->indexpos(p.k, p.l));
  p_node = node;
  node++;

  // could it be a replacement (forests are aligned)
  if (p.j > 0 && p.l > 0) {
    unsigned int t_c = 0;
    unsigned int t_rb = 0;
    getConnectedTablesRep(t, t_c, t_rb);
		R ovr = getMtrxVal(t_rb, this->f1_->over(p.i,p.j),this->f2_->over(p.k,p.l));
    // check for basepair replacement only if Algebra is of type RNA_Algebra
    if (rnaAlg_ && this->f1_->down(p.i) && this->f2_->down(p.k)) {
				R dwn = getMtrxVal(t_c, this->f1_->mdown(p.i),this->f2_->mdown(p.k));
        h_score = rnaAlg_->replacepair(this->f1_->label(p.i+1), this->f2_->label(p.k+1), dwn,
									this->f1_->label(this->f1_->getRightmostBrotherIndex2(p.i+1)),
                  this->f2_->label(this->f2_->getRightmostBrotherIndex2(p.k+1)), ovr); 
				//std::cout << "replacepair possible" << std::endl;
				if (score == h_score)
					return backtraceReplacepair(f, p, node, p_node, t_c, t_rb);
    } 
		else {
				R dwn = getMtrxVal(t_c, this->f1_->down(p.i), this->f2_->down(p.k));
				h_score = alg_->replace(this->f1_->label(p.i), dwn, this->f2_->label(p.k), ovr);
				//std::cout << "replacement possible" << std::endl;
				if (score == h_score)
					return backtraceRep(f, p, node, p_node, t_c, t_rb);
    }
	}

  // could it be a deletion
  if (p.j > 0) {
    h = p.k; // h is the node where the suffix of the split begins
    for (r = 0; r <= p.l; r++) {// for all splits of f2
      bool open = false;
      unsigned int t_c = 0;
      unsigned int t_rb = 0;
      getConnectedTablesDel(t, t_c, t_rb, open);
			L lbl = this->f1_->label(p.i);
			R dwn = getMtrxVal(t_c, this->f1_->down(p.i), this->f2_->indexpos(p.k, r));
			R ovr = getMtrxVal(t_rb, this->f1_->over(p.i, p.j), this->f2_->indexpos(h, p.l - r));
      if (open) 
        h_score =  alg_->delO(lbl, dwn, ovr);
      else 
        h_score = alg_->del(lbl, dwn, ovr);

      if (score == h_score)
        return backtraceDel(f, p, node, p_node, t, h, r);

      h = this->f2_->rb(h);
    }
  }
  // could it be an insertion
  if (p.l > 0) {
    h = p.i;
    // for all splits of f1
    for (r = 0; r <= p.j; r++) {
      bool open = false;
      unsigned int t_c = 0;
      unsigned int t_rb = 0;
      getConnectedTablesIns(t, t_c, t_rb, open);
			L lbl = this->f2_->label(p.k);
			R dwn = getMtrxVal(t_c,this->f1_->indexpos(p.i, r), this->f2_->down(p.k));
			R ovr = getMtrxVal(t_rb,this->f1_->indexpos(h, p.j - r), this->f2_->over(p.k, p.l));
      if (open) 
        h_score = alg_->insertO(dwn, lbl, ovr);
      else 
        h_score = alg_->insert(dwn, lbl, ovr);

      if (score == h_score)
        return backtraceIns(f, p, node, p_node, t, h, r);

      h = this->f1_->rb(h);
			}
  }

  std::cerr << "Strange things happening in backtrace "  <<  table_name[t] << "." << std::endl;
  std::cerr << f << "." << std::endl;
  std::exit(EXIT_FAILURE);
}

template<class R, class L, class AL>
unsigned int AlignmentAffine<R, L, AL>::backtraceReplacepair(ForestAli<L, AL> &f, CSFPair p,
    unsigned int &node, unsigned int &p_node, unsigned int t_c, unsigned int t_rb) {

      // it is a replacement
			int anchor = 0;
			if (this->anchored_) {
				if (this->f1_->getAnchor(p.i)==this->f2_->getAnchor(p.k)) {
					anchor = this->f1_->getAnchor(p.i);
					//if (anchor != 0) {
					//	std::cout << "found an anchor" << std::endl;
					//	exit(0);
					//}
				}
				else {
					std::cout << "apparently matchend unequal anchors." << std::endl;
					exit(EXIT_FAILURE);
				}
			}
			f.makeRepLabel(p_node, this->f1_->label(p.i), this->f2_->label(p.k), anchor);
      //std::cout << "it is backtrack replacepair" << std::endl;
			unsigned int rb_node;

      // set labels of leftmost child
      f.makeRepLabel(p_node+1, this->f1_->label(p.i+1), this->f2_->label(p.k+1),0);
      f.setRightBrotherIndex(p_node+1, p_node+1+1);
      f.setNumChildren(p_node+1, 0);   // base node has no children
      node++;

      // down alignment
      assert(this->f1_->noc(p.i)>=2);
      assert(this->f2_->noc(p.k)>=2);

			CSFPair down = CSFPair(p.i + 1 + 1, this->f1_->noc(p.i) - 2, p.k + 1 + 1, this->f2_->noc(p.k) - 2);
      unsigned int noc = backtrace(f, down, node, t_c);
      f.setNumChildren(p_node,noc+2);

      if (noc==0) {
          f.setRightBrotherIndex(p_node+1, p_node+1+1);
          f.setRightBrotherIndex(p_node+1+1, 0);
      }
      else {
          f.setRightBrotherIndex(f.getRightmostBrotherIndex2(p_node+1+1), node);
      }

      // set labels of rightmost child
      f.makeRepLabel(node, this->f1_->label(this->f1_->getRightmostBrotherIndex2(p.i+1+1)),
                           this->f2_->label(this->f2_->getRightmostBrotherIndex2(p.k+1+1)), 0);
      f.setRightBrotherIndex(node, 0);
      f.setNumChildren(node, 0);       // base node has no children
      node++;
      rb_node = node;                  // !!

      // right alignment
			CSFPair right = CSFPair(this->f1_->rb(p.i), p.j - 1, this->f2_->rb(p.k), p.l - 1);
      unsigned int rbs = backtrace(f, right, node, t_rb);
      if (rbs)
        f.setRightBrotherIndex(p_node, rb_node);
      else
        f.setRightBrotherIndex(p_node, 0);

      return rbs + 1;

}

template<class R, class L, class AL>
unsigned int AlignmentAffine<R, L, AL>::backtraceRep(ForestAli<L, AL> &f, CSFPair p,
    unsigned int &node, unsigned int &p_node, unsigned int t_c, unsigned int t_rb) {

			int anchor = 0;
			if (this->anchored_) {
				if (this->f1_->getAnchor(p.i)==this->f2_->getAnchor(p.k)) {
					anchor = this->f1_->getAnchor(p.i);
					//if (anchor != 0) {
					//	std::cout << "found an anchor" << std::endl;
					//	exit(0);
					//}
				}
				else {
					std::cout << "apparently matchend unequal anchors." << std::endl;
					exit(EXIT_FAILURE);
				}
			}

      // it is a replacement
      f.makeRepLabel(p_node, this->f1_->label(p.i), this->f2_->label(p.k), anchor);
      // down alignment
			CSFPair down = CSFPair(p.i + 1, this->f1_->noc(p.i), p.k + 1, this->f2_->noc(p.k));
      unsigned int noc = backtrace(f, down, node, t_c);
      f.setNumChildren(p_node, noc);
      unsigned int rb_node = node;
      // right alignment
			CSFPair right = CSFPair(this->f1_->rb(p.i), p.j - 1, this->f2_->rb(p.k), p.l - 1);
      unsigned int rbs = backtrace(f, right, node, t_rb);
      if (rbs)
        f.setRightBrotherIndex(p_node, rb_node);
      else
        f.setRightBrotherIndex(p_node, 0);

      return rbs + 1;
}

template<class R, class L, class AL>
unsigned int AlignmentAffine<R, L, AL>::backtraceDel(ForestAli<L, AL> &f,  CSFPair p,
    unsigned int &node, unsigned int &p_node,
    unsigned int t, unsigned int h, unsigned int r) {

        //std::cout << "backtraceDel from class AlignmentAffine" << std::endl;
        bool open;
        unsigned int t_c, t_rb;
        this->getConnectedTablesDel(t, t_c, t_rb, open);
        // it is a deletion
        f.makeDelLabel(p_node, this->f1_->label(p.i));
        // down alignment
				CSFPair down = CSFPair(p.i + 1, this->f1_->noc(p.i), p.k, r);
        unsigned int noc = backtrace(f, down, node, t_c);
        f.setNumChildren(p_node, noc);
        unsigned int rb_node = node;
        // right alignment
				CSFPair right = CSFPair(this->f1_->rb(p.i), p.j - 1, h, p.l - r);
        unsigned int rbs = backtrace(f, right, node, t_rb);
        if (rbs)
          f.setRightBrotherIndex(p_node, rb_node);
        else
          f.setRightBrotherIndex(p_node, 0);

        return rbs + 1;

}

template<class R, class L, class AL>
unsigned int AlignmentAffine<R, L, AL>::backtraceIns(ForestAli<L, AL> &f,  CSFPair p,
    unsigned int &node, unsigned int &p_node,
    unsigned int t, unsigned int h, unsigned int r) {

        //std::cout << "backtraceIns from class AlignmentAffine" << std::endl;
        bool open;
        unsigned int t_c, t_rb;
        this->getConnectedTablesIns(t, t_c, t_rb, open);
        // it is an insertion
        f.makeInsLabel(p_node, this->f2_->label(p.k));
        // down alignment
				CSFPair down = CSFPair(p.i, r, p.k + 1, this->f2_->noc(p.k));
        unsigned int noc = backtrace(f, down, node, t_c);
        f.setNumChildren(p_node, noc);
        unsigned int rb_node = node;
        // right alignment
				CSFPair right = CSFPair(h, p.j - r, this->f2_->rb(p.k), p.l - 1);
        unsigned int rbs = backtrace(f, right, node, t_rb);
        if (rbs)
          f.setRightBrotherIndex(p_node, rb_node);
        else
          f.setRightBrotherIndex(p_node, 0);

        return rbs + 1;
}


#endif
