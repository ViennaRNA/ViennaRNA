#ifndef _FORESTALI_H_
#define _FORESTALI_H_

#include "forest.h"

/** ForestAli is the Base class of forest Alignemnts.
 *  The pure virtual functions must be implemented when
 *  inherited to allow the construction of alignments.
 */

template <class L,class AL>
class ForestAli : public Forest<AL> {
public:
    typedef typename Forest<AL>::size_type size_type;

public:

    virtual void makeRepLabel(size_type node,L a, L b, int anchor) = 0;
    virtual void makeDelLabel(size_type node,L a) = 0;
    virtual void makeInsLabel(size_type node,L b) = 0;

public:
    ForestAli() : Forest<AL>() {};
    ForestAli(size_type size) : Forest<AL>(size) {};
};

#endif
