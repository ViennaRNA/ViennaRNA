#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include "alignment_tables.h"
#include "forestali.h"
#include <vector>
#include <sstream>

struct CSFPair {
    unsigned int i;
    unsigned int j;
    unsigned int k;
    unsigned int l;

    CSFPair(unsigned int i, unsigned int j, unsigned int k, unsigned int l) : i(i),j(j),k(k),l(l) {};
};

#define SPEEDUP true

template<class R,class L,class AL>
class Alignment {

protected:
    const Forest<L> *f1_;
    const Forest<L> *f2_;

		const bool topdown_;
		const bool anchored_;
		const bool pIndels_;
		const bool printBacktrace_;
    R localOptimum_;
    R localSubOptimum_;
    std::vector<CSFPair> localAlis_;	  // alignments already produced by getOptLocalAlignment
    double suboptimalsPercent_;         // value between 0 and 1

		void fillMatricesTopDown(bool local, bool RNA, bool speedup);
		void fillMatricesBottomUp(bool local, bool RNA, bool speedup);
		// the matrix is filled like this, implemented in affine and linear
    virtual void makeFirstCell() = 0;
    virtual void makeFirstRow() = 0;
    virtual void makeFirstCol() = 0;

		// main loops are shared by affine and linear
    void makeInnersLocal(bool RNA, bool speedup);
    void makeInnersGlobal(bool RNA, bool speedup);

		// differences in cell filling, implemented in affine and linear
    virtual void compareCSFPair(CSFPair p, bool RNA, bool speedup) = 0;
		virtual void compareCSFPairTD(CSFPair p, bool speedup) = 0;
		void recursiveFill(CSFPair startpair, bool RNA, bool speedup);
		virtual void recursiveFillAnchored(CSFPair startpair, bool RNA, bool speedup) { std::cout << "called in superclass" << std::endl;};

    virtual unsigned int backtrace(ForestAli<L,AL> &f, CSFPair p, unsigned int &node, int t=-1) = 0;

public:
		Alignment(const Forest<L> *f1, const Forest<L> *f2, const bool topdown, const bool anchored, const bool printBacktrace);
    ~Alignment();

		friend std::ostream& operator<<(std::ostream &out, const Alignment<R,L,AL> &ali) {
			ali.print(out);
			return out;
		}
		virtual void print(std::ostream &out) const = 0;

    inline void printToFile(std::string filename) const {
	    std::ofstream file (filename.c_str());
	    if (file.is_open()) {
	        print(file);
	        file.close();
	    } else {
	        std::cerr << "Unable to open file";
	        exit(EXIT_FAILURE);
	    }
		};

    inline R getLocalOptimum() {
        return localSubOptimum_;
    };

    void resetOptLocalAlignment(int suboptimalsPercent=100);

    virtual double getGlobalOptimumRelative() const = 0;
    virtual R getSILOptimum() const = 0;
    virtual void getOptGlobalAlignment(ForestAli<L,AL> &fali) = 0;
    virtual R getGlobalOptimum() const = 0;

    virtual bool nextLocalSuboptimum() = 0;
    virtual void getOptLocalAlignment(ForestAli<L,AL> &fali,unsigned int &xbasepos, unsigned int &ybasepos) = 0;
    virtual void getOptSILAlignment(ForestAli<L,AL> &fali,unsigned int &ybasepos) = 0;

		virtual bool computed(const unsigned long i, const unsigned long j) const = 0; 
		virtual void setComputed(const unsigned long i, const unsigned long j) = 0; 

		// Index transition functions to subproblems
	
    CSFPair downReplacepair(CSFPair p) const {
      return CSFPair(p.i + 1 + 1, this->f1_->noc(p.i) - 2, p.k + 1 + 1, this->f2_->noc(p.k) - 2);  // !!    mdown !!
    }


    CSFPair downDeletePairOnly(CSFPair p, unsigned int r) const {
      assert (r-2 >= 0);
      return CSFPair(p.i + 1 + 1, this->f1_->noc(p.i) - 2, p.k + 1, r - 2);
    }
    CSFPair downDeletePairAndBases(CSFPair p, unsigned int r) const {
      return CSFPair(p.i + 1 + 1, this->f1_->noc(p.i) - 2, p.k, r);
    }
    CSFPair downDeletePairAndLeftBase(CSFPair p, unsigned int r) const {
      assert (r-1 >= 0);
      return CSFPair(p.i + 1 + 1, this->f1_->noc(p.i) - 2, p.k, r - 1);
    }
    CSFPair downDeletePairAndRightBase(CSFPair p, unsigned int r) const {
      assert (r-1 >= 0);
      return CSFPair(p.i + 1 + 1, this->f1_->noc(p.i) - 2, p.k + 1, r-1);
    }


    CSFPair downInsertPairOnly(CSFPair p, unsigned int r) const {
      assert (r-2 >= 0);
      return CSFPair(p.i + 1, r - 2, p.k + 1 + 1, this->f2_->noc(p.k) - 2);
    }
    CSFPair downInsertPairAndBases(CSFPair p, unsigned int r) const {
      return CSFPair(p.i, r, p.k + 1 + 1, this->f2_->noc(p.k) - 2);
    }
    CSFPair downInsertPairAndLeftBase(CSFPair p, unsigned int r) const {
      assert (r-1 >= 0);
      return CSFPair(p.i, r - 1, p.k + 1 + 1, this->f2_->noc(p.k) - 2);
    }
    CSFPair downInsertPairAndRightBase(CSFPair p, unsigned int r) const {
      assert (r-1 >= 0);
      return CSFPair(p.i + 1, r - 1, p.k + 1 + 1, this->f2_->noc(p.k) - 2);
    }


    CSFPair downReplace(CSFPair p) const {
      return CSFPair(p.i + 1, this->f1_->noc(p.i), p.k + 1, this->f2_->noc(p.k));
    }
    CSFPair rightReplace(CSFPair p) const {
      return CSFPair(this->f1_->rb(p.i), p.j - 1, this->f2_->rb(p.k), p.l - 1);
    }

    CSFPair downDelete(CSFPair p, unsigned int r, unsigned int h) const {
      return CSFPair(p.i + 1, this->f1_->noc(p.i), p.k, r);
    }
    CSFPair rightDelete(CSFPair p, unsigned int r, unsigned int h) const {
      return CSFPair(this->f1_->rb(p.i), p.j - 1, h, p.l - r);
    }

    CSFPair downInsert(CSFPair p, unsigned int r, unsigned int h) const {
      return CSFPair(p.i, r, p.k + 1, this->f2_->noc(p.k));
    }
    CSFPair rightInsert(CSFPair p, unsigned int r, unsigned int h) const {
      return CSFPair(h, p.j - r, this->f2_->rb(p.k), p.l - 1);
    }



};

template<class R,class L,class AL>
class AlignmentLinear : public Alignment<R,L,AL> {
	private:
		TAD_DP_TableLinear<R> *mtrx_;
    const Algebra<R,L> *alg_;
    const RNA_Algebra<R,L> *rnaAlg_;

    void compareCSFPair(CSFPair p, bool RNA, bool speedup);
		void compareCSFPairTD(CSFPair p, bool speedup);

	public:

		void print(std::ostream &out) const { out << "linear ali's matrix" << std::endl << *mtrx_; };

    AlignmentLinear(const Forest<L> *f1, const Forest<L> *f2,const Algebra<R,L> &alg, const bool topdown, const bool anchored, bool local, bool printBacktrace, bool speedup=SPEEDUP);
    AlignmentLinear(const Forest<L> *f1, const Forest<L> *f2,const RNA_Algebra<R,L> &rnaAlg, const bool topdown, const bool anchored, bool local, bool printBacktrace, bool speedup=SPEEDUP);
    void makeFirstCell();
    void makeFirstRow();
    void makeFirstCol();
		void recursiveFillAnchored(CSFPair startpair, bool RNA, bool speedup);

    unsigned int backtrace(ForestAli<L,AL> &f, CSFPair p, unsigned int &node, int t=-1);

    inline R getMtrxVal(unsigned long i,unsigned long j) const {
			return mtrx_->getMtrxVal(i,j);
    }

		// TODO hier wird gar nicht darauf geachted ob alg oder rnaalg
    inline void setMtrxVal(unsigned long i,unsigned long j,R val) {
			mtrx_->setMtrxVal(i,j,val);
      // here we calculate local similarity on the fly
      this->localOptimum_ = alg_->choice(this->localOptimum_, val);
    }

    // function should only be used with similarity based algebras
    double getGlobalOptimumRelative() const;
    R getSILOptimum() const;
    void getOptGlobalAlignment(ForestAli<L,AL> &fali);
    R getGlobalOptimum() const;

    bool nextLocalSuboptimum();
    void getOptLocalAlignment(ForestAli<L,AL> &fali,unsigned int &xbasepos, unsigned int &ybasepos);
    void getOptSILAlignment(ForestAli<L,AL> &fali,unsigned int &ybasepos);

		bool computed(const unsigned long i, const unsigned long j) const { return mtrx_->computed(i,j); }; 
		void setComputed(const unsigned long i, const unsigned long j) { mtrx_->setComputed(i,j); }; 

    // virtual, for replacepair
    virtual inline R computeReplacementScore(CSFPair p, std::string & backtrack_as) const {
        backtrack_as = 'R';
        return alg_->replace(this->f1_->label(p.i),
          getMtrxVal(this->f1_->down(p.i), this->f2_->down(p.k)),
          this->f2_->label(p.k),
          getMtrxVal(this->f1_->over(p.i, p.j), this->f2_->over(p.k, p.l)));
    }

    // for deletepair
    virtual inline R computeDeletionScore(CSFPair p, std::string & backtrack_as, unsigned int r, unsigned int h) const {
      backtrack_as = 'D';
      return alg_->del(this->f1_->label(p.i),
          getMtrxVal(this->f1_->down(p.i), this->f2_->indexpos(p.k, r)),
          getMtrxVal(this->f1_->over(p.i, p.j), this->f2_->indexpos(h, p.l - r)));
    }

    // for insertpair
    virtual inline R computeInsertionScore(CSFPair p, std::string & backtrack_as, unsigned int r, unsigned    int h) const {
      backtrack_as = 'I';
      return alg_->insert(getMtrxVal(this->f1_->indexpos(p.i, r), this->f2_->down(p.k)),
        this->f2_->label(p.k),
        getMtrxVal(this->f1_->indexpos(h, p.j - r), this->f2_->over(p.k, p.l)));
    }

    // for replacepair
    virtual inline R computeReplacementScoreRNA(CSFPair p, std::string & backtrack_as) const {
      if ((this->f1_->isInternalNode(p.i)) 
					&& (this->f2_->isInternalNode(p.k))) {
        backtrack_as = 'P';
        return rnaAlg_->replacepair(this->f1_->label(p.i+1),
          this->f2_->label(p.k+1),
          getMtrxVal(this->f1_->mdown(p.i),this->f2_->mdown(p.k)),
          this->f1_->label(this->f1_->getRightmostBrotherIndex(p.i+1)),
          this->f2_->label(this->f2_->getRightmostBrotherIndex(p.k+1)),
          getMtrxVal(this->f1_->over(p.i,p.j),this->f2_->over(p.k,p.l)));
      }
      else if ((this->f1_->isInternalNode(p.i))) {
        backtrack_as = 'x';
        return rnaAlg_->worst_score();
      }
      else if ((this->f2_->isInternalNode(p.k))) {
        backtrack_as = 'x';
        return rnaAlg_->worst_score();
      }
      else {
        backtrack_as = 'R';
        return rnaAlg_->replace(this->f1_->label(p.i),
          getMtrxVal(this->f1_->down(p.i), this->f2_->down(p.k)),
          this->f2_->label(p.k),
          getMtrxVal(this->f1_->over(p.i, p.j), this->f2_->over(p.k, p.l)));
      }
    }

    // for deletepair
    virtual inline R computeDeletionScoreRNA(CSFPair p, std::string & backtrack_as, unsigned int r, unsigned   int h) const {
      backtrack_as = 'D';
      if ((this->f1_->isInternalNode(p.i))) {
        
        R score = rnaAlg_->worst_score();
        
        //std::cout << "computing score for subforest " << p.i << ", " << p.j << ", " << p.k << ", " << p.l << std::endl; 
        //std::cout << "r = " << r << std::endl;
        //std::cout << "exit in computeDeletionScore" << std::endl;
        //exit(0);
    
        unsigned int la = 0, ra = 0, lb = 0, rb = 0;
        lb = p.k;

        // alignment of children:
        // delete both
        //assert enough leaves - we need none
        if ((r==0 || r==1 || r>1)) {
         R h_score = compute_delete_pairing_and_bases(p, r, h);
         score = rnaAlg_->choice(score, h_score);
         //std::cout << "computed score delete_pairing_and_bases " << h_score << std::endl;
        }
        // delete left
        //assert enough leaves - one to the right
        if (r==1 || r>1) {
          rb = this->f2_->rb(p.k, r-1);
          if  (this->f2_->isLeaf(rb)) {
            R h_score = compute_delete_pairing_and_left_base(p, r, h);
            score = rnaAlg_->choice(score, h_score);
            //std::cout << "computed score delete_pairing_and_left_base " << h_score << std::endl;
          }
        }
        // delete right 
        //assert enough leaves - one to the left
        if ((r==1 || r>1) && (this->f2_->isLeaf(lb))) {
         R h_score = compute_delete_pairing_and_right_base(p, r, h);
         score = rnaAlg_->choice(score, h_score);
         //std::cout << "computed score delete_pairing_and_right_base " << h_score << std::endl;
        }
        // delete none - replace both
        //assert enough leaves - one to the left and one to the right
        if (r>1) { 
          rb = this->f2_->rb(p.k, r-1);
          if (this->f2_->isLeaf(lb) && this->f2_->isLeaf(rb)) {
            R h_score = compute_delete_pairing_only(p, r, h);
            score = rnaAlg_->choice(score, h_score);
            //std::cout << "computed score delete_pairing_only " << h_score << std::endl;
          }
        }
        return score;
      }
      //std::cout << "computed score delete " << rnaAlg_->del(this->f1_->label(p.i),
      //    getMtrxVal(this->f1_->down(p.i), this->f2_->indexpos(p.k, r)),
      //    getMtrxVal(this->f1_->over(p.i, p.j), this->f2_->indexpos(h, p.l - r))) << std::endl;
      return rnaAlg_->del(this->f1_->label(p.i),
          getMtrxVal(this->f1_->down(p.i), this->f2_->indexpos(p.k, r)),
          getMtrxVal(this->f1_->over(p.i, p.j), this->f2_->indexpos(h, p.l - r)));
    }

    // for insertpair
    virtual inline R computeInsertionScoreRNA(CSFPair p, std::string & backtrack_as, unsigned int r, unsigned  int h) const {
      std::ostringstream backtrack_as_oss;
      backtrack_as_oss << "I";
      if ((this->f2_->isInternalNode(p.k))) {
        
        R score = rnaAlg_->worst_score();
        
        unsigned int la = 0, ra = 0, lb = 0, rb = 0;
        la = p.i;

        // alignment of children:
        // insert both
        //assert enough leaves - we need none
        if ((r==0 || r==1 || r>1)) {
          R h_score = compute_insert_pairing_and_bases(p, r, h);
          score = rnaAlg_->choice(score, h_score);
          backtrack_as_oss << "Ip&b " << score << " ";
        }
        // insert left
        //assert enough leaves - one to the right
        if (r==1 || r>1) {
          ra = this->f1_->rb(p.i, r-1);
          if (this->f1_->isLeaf(ra)) {
            R h_score = compute_insert_pairing_and_left_base(p, r, h);
            score = rnaAlg_->choice(score, h_score);
            if (score == h_score)
              backtrack_as_oss << "Ip&l " << score << " ";
          }
        }
        // insert none - replace both
        //assert enough leaves - one to the left and one to the right
        if (r>1) {
          ra = this->f1_->rb(p.i, r-1);
          if (this->f1_->isLeaf(la) && this->f1_->isLeaf(ra)) {
            R h_score = compute_insert_pairing_only(p, r, h);
            score = rnaAlg_->choice(score, h_score);
            //std::cout << "computed score insert_pairing_only " << h_score << std::endl;
            if (score == h_score)
              backtrack_as_oss << "IpO " << score << " ";
          }
        }
        // insert right 
        //assert enough leaves - one to the left
        if (r==1 || r>1) {
          ra = this->f1_->rb(p.i, r-1);
          //std::cout << "right in a " << ra << std::endl;
          if (this->f1_->isLeaf(la) && this->f1_->isLeaf(ra)) {
            R h_score = compute_insert_pairing_and_right_base(p, r, h);
            score = rnaAlg_->choice(score, h_score);
            //std::cout << "computed score insert_pairing_and_right_base " << h_score << std::endl;
            if (score == h_score)
              backtrack_as_oss << "Ip&R " << score << " ";
          }
        }

        backtrack_as = backtrack_as_oss.str();
        return score;
      }
      //std::cout << "computed score insert " << rnaAlg_->insert(getMtrxVal(this->f1_->indexpos(p.i, r),     this->f2_->down(p.k)),
		}

		// Pair deletion functions
    R compute_delete_pairing_only(CSFPair p, unsigned int r, unsigned int h) const {
      CSFPair mdown_p = this->downDeletePairOnly(p, r);
      R mdown = getMtrxVal(this->f1_->indexpos(mdown_p.i, mdown_p.j), this->f2_->indexpos(mdown_p.k, mdown_p. l));
      R over = getMtrxVal(this->f1_->over(p.i, p.j), this->f2_->indexpos(h, p.l - r));
      return rnaAlg_->deletePairOnly(this->f1_->label(p.i+1), this->f2_->label(p.k), 
          mdown,
          this->f1_->label(this->f1_->getRightmostBrotherIndex(p.i+1)),
          this->f2_->label(this->f2_->rb(p.k, r)),
          over);
    }
  
    R compute_delete_pairing_and_right_base(CSFPair p, unsigned int r, unsigned int h) const {
      CSFPair mdown_p = this->downDeletePairAndRightBase(p, r);
      R mdown = getMtrxVal(this->f1_->indexpos(mdown_p.i, mdown_p.j), this->f2_->indexpos(mdown_p.k, mdown_p. l));
      R over = getMtrxVal(this->f1_->over(p.i, p.j), this->f2_->indexpos(h, p.l - r));
      return rnaAlg_->deletePairAndRightBase(this->f1_->label(p.i+1), this->f2_->label(p.k),
        mdown,
        this->f1_->label(this->f1_->getRightmostBrotherIndex(p.i+1)),
        over);
    }
  
    R compute_delete_pairing_and_left_base(CSFPair p, unsigned int r, unsigned int h) const {
      CSFPair mdown_p = this->downDeletePairAndLeftBase(p, r);
      R mdown = getMtrxVal(this->f1_->indexpos(mdown_p.i, mdown_p.j), this->f2_->indexpos(mdown_p.k, mdown_p. l));
      R over = getMtrxVal(this->f1_->over(p.i, p.j), this->f2_->indexpos(h, p.l - r));
      return rnaAlg_->deletePairAndLeftBase(this->f1_->label(p.i+1),
          mdown,
          this->f1_->label(this->f1_->getRightmostBrotherIndex(p.i+1)),
          this->f2_->label(this->f2_->rb(p.k, r)),
          over);
    }
  
    R compute_delete_pairing_and_bases(CSFPair p, unsigned int r, unsigned int h) const {
      CSFPair mdown_p = this->downDeletePairAndBases(p, r);
      R mdown = getMtrxVal(this->f1_->indexpos(mdown_p.i, mdown_p.j), this->f2_->indexpos(mdown_p.k, mdown_p. l));
      R over = getMtrxVal(this->f1_->over(p.i, p.j), this->f2_->indexpos(h, p.l - r));
      return rnaAlg_->deletePairAndBases(this->f1_->label(p.i+1),
        mdown,
        this->f1_->label(this->f1_->getRightmostBrotherIndex(p.i+1)),
        over);
    }

		// Pair insertion functions
    R compute_insert_pairing_only(CSFPair p, unsigned int r, unsigned int h) const {
       CSFPair mdown_p = this->downInsertPairOnly(p, r);
       R mdown = getMtrxVal(this->f1_->indexpos(mdown_p.i, mdown_p.j), this->f2_->indexpos(mdown_p.k, mdown_p.l));
       R over = getMtrxVal(this->f1_->indexpos(h, p.j - r), this->f2_->over(p.k, p.l));
       return rnaAlg_->insertPairOnly(this->f1_->label(p.i), this->f2_->label(p.k+1), 
         mdown,
         this->f1_->label(this->f1_->rb(p.i, r)),
         this->f2_->label(this->f2_->getRightmostBrotherIndex(p.k+1)),
         over);
    }

    R compute_insert_pairing_and_right_base(CSFPair p, unsigned int r, unsigned int h) const {
       CSFPair mdown_p = this->downInsertPairAndRightBase(p, r);
       R mdown = getMtrxVal(this->f1_->indexpos(mdown_p.i, mdown_p.j), this->f2_->indexpos(mdown_p.k, mdown_p.l));
       R over = getMtrxVal(this->f1_->indexpos(h, p.j - r), this->f2_->over(p.k, p.l));
       return rnaAlg_->insertPairAndRightBase(this->f1_->label(p.i), this->f2_->label(p.k+1),
         mdown,
         this->f2_->label(this->f2_->getRightmostBrotherIndex(p.k+1)),
         over);
    }

    R compute_insert_pairing_and_left_base(CSFPair p, unsigned int r, unsigned int h) const {
       CSFPair mdown_p = this->downInsertPairAndLeftBase(p, r);
       R mdown = getMtrxVal(this->f1_->indexpos(mdown_p.i, mdown_p.j), this->f2_->indexpos(mdown_p.k, mdown_p.l));
       R over = getMtrxVal(this->f1_->indexpos(h, p.j - r), this->f2_->over(p.k, p.l));
       return rnaAlg_->insertPairAndLeftBase(this->f2_->label(p.k+1),
         mdown,
         this->f1_->label(this->f1_->rb(p.i, r)),
         this->f2_->label(this->f2_->getRightmostBrotherIndex(p.k+1)),
         over);
    }

    R compute_insert_pairing_and_bases(CSFPair p, unsigned int r, unsigned int h) const {
        CSFPair mdown_p = this->downInsertPairAndBases(p, r);
        R mdown = getMtrxVal(this->f1_->indexpos(mdown_p.i, mdown_p.j), this->f2_->indexpos(mdown_p.k,        mdown_p.l));
        R over = getMtrxVal(this->f1_->indexpos(h, p.j - r), this->f2_->over(p.k, p.l));
        return rnaAlg_->insertPairAndBases(this->f2_->label(p.k+1),
          mdown,
          this->f2_->label(this->f2_->getRightmostBrotherIndex(p.k+1)),
          over);
    }

};

template<class R,class L,class AL>
class AlignmentAffine : public Alignment<R,L,AL> {
	private:
		TAD_DP_TableAffine<R> *mtrx_;
    const AlgebraAffine<R,L> *alg_;
    const RNA_AlgebraAffine<R,L> *rnaAlg_;
		int localOptimumTable_;

    void compareCSFPair(CSFPair p, bool RNA, bool speedup);
		void compareCSFPairTD(CSFPair p, bool speedup);

		// special for the seven tables - maybe put in the tables class or outside class?
		void getConnectedTablesRep(const unsigned int t, unsigned int &t_c, unsigned int &t_rb);
		void getConnectedTablesIns(const unsigned int t, unsigned int &t_c, unsigned int &t_rb, bool &open);
		void getConnectedTablesDel(const unsigned int t, unsigned int &t_c, unsigned int &t_rb, bool &open);

    unsigned int backtraceReplacepair(ForestAli<L, AL> &f, CSFPair p, unsigned int &node, unsigned int &p_node, const unsigned int t_c, const unsigned int t_rb);
    unsigned int backtraceRep(ForestAli<L, AL> &f, CSFPair p, unsigned int &node, unsigned int &p_node, const unsigned int t_c, const unsigned int t_rb);
    unsigned int backtraceIns(ForestAli<L, AL> &f, CSFPair p, unsigned int &node, unsigned int &p_node, const unsigned int t, unsigned int h, unsigned int r);
    unsigned int backtraceDel(ForestAli<L, AL> &f, CSFPair p, unsigned int &node, unsigned int &p_node, const unsigned int t, unsigned int h, unsigned int r);
	public:

		void print(std::ostream &out) const { out << "affine ali's matrix" << std::endl << *mtrx_; };

    AlignmentAffine(const Forest<L> *f1, const Forest<L> *f2,const AlgebraAffine<R,L> &alg, 
				const bool topdown, const bool anchored, bool local, bool printBacktrace, bool speedup=SPEEDUP);
    AlignmentAffine(const Forest<L> *f1, const Forest<L> *f2,const RNA_AlgebraAffine<R,L> &rnaAlg, 
				const bool topdown, const bool anchored, bool local, bool printBacktrace, bool speedup=SPEEDUP);
    void makeFirstCell();
    void makeFirstRow();
    void makeFirstCol();

		void recursiveFillAnchored(CSFPair startpair, bool RNA, bool speedup);

    inline R getMtrxVal(const int t, const unsigned long i, const unsigned long j) const {
			return mtrx_->getMtrxVal(t,i,j);
    }
    inline void setMtrxVal(const int t, const unsigned long i, const unsigned long j, R val) {
			mtrx_->setMtrxVal(t,i,j,val);
			// NOTE: we use alg_ even for RNA, but as we use the choice function only, that's ok
      // here we calculate local similarity on the fly
      this->localOptimum_ = alg_->choice(this->localOptimum_, val);
      if (this->localOptimum_ == val)
        	localOptimumTable_ = t;
    }

    // function should only be used with similarity based algebras
    double getGlobalOptimumRelative() const;
    R getSILOptimum() const;
    void getOptGlobalAlignment(ForestAli<L,AL> &fali);
    R getGlobalOptimum() const;

    bool nextLocalSuboptimum();
    void getOptLocalAlignment(ForestAli<L,AL> &fali,unsigned int &xbasepos, unsigned int &ybasepos);
    void getOptSILAlignment(ForestAli<L,AL> &fali,unsigned int &ybasepos);
    void foundOptLocalAlignment(ForestAli<L, AL> &fali, CSFPair csfp, unsigned int &start1, unsigned int &start2, unsigned int &end1, unsigned int &end2);
    unsigned int backtrace(ForestAli<L,AL> &f, CSFPair p, unsigned int &node, int t=-1);

		bool computed(const unsigned long i, const unsigned long j) const { return mtrx_->computed(i,j); }; 
		void setComputed(const unsigned long i, const unsigned long j) { mtrx_->setComputed(i,j); }; 
};

#endif
