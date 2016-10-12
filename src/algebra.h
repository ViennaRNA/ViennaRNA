#ifndef _ALGEBRA_H_
#define _ALGEBRA_H_

// Algebra is a virtual class (interface) for algebra used by the Alignment template classes.
// This is in the spirit of the Algebraic Dynamic Programming (ADP) approach where a Dynamic Programming
// algorithm is seperated into a grammar and an algebra.

template<class R, class L>
class Algebra {
public:
    virtual R empty() const =0;								/**< Result for the empty tree alignment */
    virtual R replace(L a,R down, L b, R over) const =0;	/**< Result for the tree edit function 'replace' */
    virtual R del(L a,R down, R over) const =0;				/**< Result for the tree edit function 'delete' */
    virtual R insert(R down,L b,R over) const =0;			/**< Result for the tree edit function 'insert' */
    virtual R choice(R a,R b) const =0;						/**< The choice function. Commonly used functions are 'min' and 'max' */
    virtual R worst_score() const =0;						/**< The worst_score with respect to choice is specified by this function */

    virtual ~Algebra() {};
};

template<class R, class L>
class AlgebraAffine : public Algebra<R,L> {
	public:
    virtual R delO(L a,R down, R over) const =0;				/**< Result for the tree edit function 'delete' */
    virtual R insertO(R down,L b,R over) const =0;			/**< Result for the tree edit function 'insert' */
};

// Extended Algebra for aligning RNA secondary structure trees where basepair replacements
// are considered as a single edit operation.
// -> Extra function for replacing a pair and its bases.
// Pair indel und Pair match are handled within the usual indel and match function (check for ALPHA_BASEPAIR).

template<class R, class L>
class RNA_Algebra : public Algebra<R,L> {
public:
    /** Result for the replacement of a basepair */
    virtual R replacepair(L la, L lb, R down, L ra, L rb, R over) const =0;

    virtual R deletePairOnly(L la, L lb, R down, L ra, L rb, R over) const = 0;
    virtual R deletePairAndBases(L la, R down, L ra, R over) const = 0;
    virtual R deletePairAndLeftBase(L la, R down, L ra, L rb, R over) const = 0;
    virtual R deletePairAndRightBase(L la, L lb, R down, L ra, R over) const = 0;

    virtual R insertPairOnly(L la, L lb, R down, L ra, L rb, R over) const = 0;
    virtual R insertPairAndBases(L lb, R down, L rb, R over) const = 0;
    virtual R insertPairAndLeftBase(L lb, R down, L ra, L rb, R over) const = 0;
    virtual R insertPairAndRightBase(L la, L lb, R down, L rb, R over) const = 0;
};

template<class R, class L>
class RNA_AlgebraAffine : public AlgebraAffine<R,L> {
public:
    /** Result for the replacement of a basepair */
    virtual R replacepair(L la, L lb, R down, L ra, L rb, R over) const =0;

    virtual R deletePairOnly(L la, L lb, R down, L ra, L rb, R over) const = 0;
    virtual R deletePairAndBases(L la, R down, L ra, R over) const = 0;
    virtual R deletePairAndLeftBase(L la, R down, L ra, L rb, R over) const = 0;
    virtual R deletePairAndRightBase(L la, L lb, R down, L ra, R over) const = 0;

    virtual R insertPairOnly(L la, L lb, R down, L ra, L rb, R over) const = 0;
    virtual R insertPairAndBases(L lb, R down, L rb, R over) const = 0;
    virtual R insertPairAndLeftBase(L lb, R down, L ra, L rb, R over) const = 0;
    virtual R insertPairAndRightBase(L la, L lb, R down, L rb, R over) const = 0;

};



// Algebra for the tree edit model

template<class R, class L>
class SZAlgebra {
public:
    virtual R empty() const =0;								/**< Result for the empty tree alignment */
    virtual R replace(L a,R down, L b) const =0;			/**< Result for the tree edit function 'replace' */
    virtual R del(L a,R down) const =0;						/**< Result for the tree edit function 'delete' */
    virtual R insert(R down,L b) const =0;					/**< Result for the tree edit function 'insert' */
    virtual R choice(R a,R b) const =0;						/**< The choice function. Commonly used functions are 'min' and 'max' */
    virtual R worst_score() const =0;						/**< The worst_score with respect to choice is specified by this function */

    virtual ~SZAlgebra() {};
};

#endif





