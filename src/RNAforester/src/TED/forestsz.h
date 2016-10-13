#ifndef _FORESTSZ_H
#define _FORESTSZ_H

#include <cassert>
#include <iostream>
#include "misc.h"



template <class L>
class ForestSZ {
protected:
    unsigned int size_;
    unsigned int *lml_;         // postorder index of leftmost leaf descandant of the subtree rootet at T[i]
    L *lb_;			   // labels of nodes in postorder
    bool *keyroot_;	   // is node a keyroot

    void calcKeyroots();

public:
    typedef unsigned int size_type;
    typedef L label_type;

    ForestSZ() : size_(0), lml_(NULL), lb_(NULL),keyroot_(NULL) {};
    ForestSZ(unsigned int nrOfNodes);
    ForestSZ(const ForestSZ<L> &f);
    ~ForestSZ();

    inline size_type size() const {
        return size_;
    };
    inline unsigned int lml(unsigned long i) const {
        return lml_[i];
    };
    inline L label(unsigned long i) const {
        return lb_[i];
    };
    inline bool keyroot(unsigned long i) const {
        return keyroot_[i];
    };

    unsigned int numKeyroots() const {
        unsigned int count=0;
        for (unsigned int i=0; i<size_; i++) {
            if (keyroot_[i])
                count++;
        }
        return count;
    }

    void printParameters(const std::string &name) const {
        std::cout << name.c_str() << "# size: " << size() << std::endl;
        std::cout << name.c_str() << "# keyroots: " << numKeyroots() << std::endl;
    }
};

#endif

