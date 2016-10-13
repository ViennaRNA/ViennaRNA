#ifndef _FOREST_BASE_H
#define _FOREST_BASE_H

#include <algorithm>
#include <cassert>
#include "misc.h"



/** ForestBase is the base class of the template class Forest<L>.
 *  To reduce the size of compiled programs functions and variables that are
 *  independent of the labelling are implemented in this class.
 */
class ForestBase {

public:
    typedef unsigned int size_type;

private:
    bool isSumUpCSFValid_;
    bool isRMBValid_;

    size_type getNumRightBrothers(size_type i) const;
    size_type countLeaves(size_type i,size_type k) const;

protected:
    size_type size_;         /**< size is the number of nodes of a tree */
    size_type *rb_;          /**< rb_[i] stores the preorder index of the rightbrother node of the ith node, or 0 if there is none. */
    size_type *noc_;         /**< noc_[i] stores the  number of children of the ith node */
    size_type *sumUpCSF_;    /**< sumUpCSF_[i] stores the sum of non empty closed subforests of nodes k<i. This is used to map index pairs to single indexes. */
    size_type *rmb_;         /**< rmb_[i] stores the index of the rightmost brother of the ith node. */
    size_type *anchors_;		 /**< anchors_[i] marks each node with an anchor number or 0. */

    // Alignment construction functions for backtrace routine of class Alignment

    /** Allocate memory and initialize variables */
    void initialize(size_type size, bool anchored);



public:
		// TODO war nicht immer public
    inline void setSize(size_type size) {
        size_=size;
    };

    // TODO diese anchor funktionen, sind sie noetig?
    // set the anchor of a node
    inline void setAnchor(const size_type node, const size_type index) { //assert(node <= size_); 
      anchors_[node] = index; }
    inline int getAnchor(const size_type node) const { //assert(node <= size_); 
      return anchors_[node]; }
    //inline std::vector<int> getAnchors() const { return anchors_; }

    /** Calculate sumUpCSF_ from rb_ and noc_. */
    void calcSumUpCSF();
    /** Calculate rmb_ from rb_ and noc_ */
    void calcRMB();

    inline void setNumChildren(size_type node,size_type num) {
        noc_[node]=num;
    };
    inline void setRightBrotherIndex(size_type node,size_type index) {
        rb_[node]=index;
    };
		// TODO ende

    /** Default constructor creates an empty forest */
    ForestBase() : isSumUpCSFValid_(false), isRMBValid_(false), size_(0),rb_(NULL),noc_(NULL),sumUpCSF_(NULL),rmb_(NULL) {};
    /** Allocates the memory for a forest of "size" nodes */
    ForestBase(size_type size);
    /** Copy constructor */
    ForestBase(const ForestBase &fBase);
    /** Frees allocated memory */
    ~ForestBase();

    inline size_type size() const {
        return size_;
    };               /**< Returns size_ */
    inline size_type noc(size_type i) const {
        return noc_[i];
    };   /**< Returns noc_[i] */
    inline size_type rb(size_type i) const {
        return rb_[i];
    };     /**< Returns rb_[i] */
    size_type rb(size_type i, size_type k) const;                 /**< Returns the kth brother of i */


    inline bool isLeaf(size_type i) const {
        if (noc_[i])
            return false;
        else
            return true;
    };

    inline bool isInternalNode(size_type i) const {
        return !isLeaf(i);
    };

    /** Calculates the maximal length of a closed subforest with node index i in constant time.
     *  len(i)=sumUpCSF[i+1]-sumUpCSF[i], if i<size(F) and 1 otherwise.
     *  Requires that function calcSumUpCSF() was already called.
     */
    inline size_type getMaxLength(size_type i) const {
        assert(isSumUpCSFValid_);

        if (size_==1)
            return 1;
        else
            return sumUpCSF_[i+1]-sumUpCSF_[i]; /* len(i)=sumUpCSF[i+1]-sumUpCSF[i], if i<nodes(F)*/
    }

    /** Returns the number of all non empty closed subforests in a forest
     *  Requires that function calcSumUpCSF() was already called.
     */
    inline size_type getNumCSFs() const {
        assert(isSumUpCSFValid_);
        return sumUpCSF_[size_];  // max node is size_-1, hence sumUpCSF_[size_] stores the number of all csfs
    };

    /** Returns the index if the rightmost brother node of i in constant time.
     *  Requires that function calcRMB() was already called.
     */
    inline size_type getRightmostBrotherIndex(size_type i) const {
        assert(isRMBValid_);
        assert(i<size_);
        return rmb_[i];
    };

    /** Returns the index if the rightmost brother node of i.
     *  Does not require that function calcRMB() was already called and runs linear in the number of brother nodes.
     */
    inline size_type getRightmostBrotherIndex2(size_type i) const {
        size_type h=i;
        while (rb_[h])
            h=rb_[h];

        return h;
    };

    size_type maxDegree() const;
    size_type numLeaves() const;
    size_type maxDepth() const;

    /** Returns the number of leaf nodes until the leftmost leave node of i.
     *  This is useful to calculate the nucleotide positions of local alignments.
     */
    size_type countLeaves(size_type i) const;

    /** @name Index transition functions */
    //@{

    /** Calculates one dimensional index for a closed subforest (i,j) */
    inline size_type indexpos(size_type i,size_type j) const {
        assert(isSumUpCSFValid_);

        if (j==0)
            return 0;
        else
            return sumUpCSF_[i]+j-1;
    }

    /** Returns indexpos(i+1,noc(i)) for a csf (i,j) */
    inline size_type down(size_type i) const {
        return indexpos(i+1,noc_[i]);
    };
    /** Returns indexpos(rb(i),j-1) for a csf (i,j) */
    inline size_type over(size_type i,size_type j) const {
        return indexpos(rb_[i],j-1);
    }
    /** Returns indexpos(rb(i),getMaxLength(i)-1) for a node index i */
    inline size_type over(size_type i) const {
        return indexpos(rb_[i],getMaxLength(i)-1);
    }

    /** Returns the mapped index of csf (i2,..,in-1) where i1,...,in are the children of i.
      * This transition is important for the extended alignment model for RNA structures where
      * a P-node and the pairing bases can be replaced by a single edit operation
      */
    inline size_type mdown(size_type i) const {
        if (noc_[i]<=2)
            return 0;
        else
            return indexpos(i+1+1,noc_[i]-2);
    };

    //@}


    /** Tests if (i,j) and (i2,j2) are disjoint.
     *  Intuitively, (i,j) and (i2,j2) are disjoint if they do not include
     *  each other and dont intersect.
     */
    inline bool isDisjoint(size_type i, size_type j, size_type i2, size_type j2) const {
        size_type max_node;

        // the empty forest is included in any forest
        if (j==0 || j2==0)
            return false;

        // (i2,j2) included in (i,j)

        max_node=getRightmostBrotherIndex(i);
        if (noc(max_node))
            max_node=getRightmostBrotherIndex(max_node+1);

        if (i2>=i && i2<=max_node)
            return false;

        if (rb(i2,j2-1)>=i && rb(i2,j2-1)<=max_node)
            return false;

        // (i,j) included in (i2,j2)
        max_node=getRightmostBrotherIndex(i2);
        if (noc(max_node))
            max_node=getRightmostBrotherIndex(max_node+1);


        if (i>=i2 && i<=max_node)
            return false;

        if (rb(i,j-1)>=i2 && rb(i,j-1)<=max_node)
            return false;

        return true;
    }
};

#endif
