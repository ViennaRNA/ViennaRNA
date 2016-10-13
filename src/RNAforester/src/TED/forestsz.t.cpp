#ifndef _FORESTSZ_T_CPP_
#define _FORESTSZ_T_CPP_

#include "misc.h"
#include "forestsz.h"

#include <map>
#include <string.h>

// Forest<T>

template <class L>
ForestSZ<L>::ForestSZ(unsigned int nrOfNodes)
        : size_(nrOfNodes) {
    lml_=new unsigned int[nrOfNodes];
    lb_=new L[nrOfNodes];
    keyroot_=new bool[nrOfNodes];
}

template <class L>
ForestSZ<L>::ForestSZ(const ForestSZ<L> &f) {
    size_=f.size();

    lml_=new unsigned int[f.size()];
    lb_=new L[f.size()];
    keyroot_=new bool[f.size()];

    memcpy(lml_,f.lml_,sizeof(unsigned int)*size_);
    memcpy(keyroot_,f.keyroot_,sizeof(bool)*size_);

    for (unsigned int i=0; i<size_; i++)
        lb_[i]=f.lb_[i];
}

template <class L>
ForestSZ<L>::~ForestSZ() {
    DELETE_ARRAY(lml_);
    DELETE_ARRAY(lb_);
    DELETE_ARRAY(keyroot_);
}

template <class L>
void ForestSZ<L>::calcKeyroots() {
    std::map<unsigned int,unsigned int> keyrootMap;

    for (unsigned int i=0; i<size_; i++) {
        keyroot_[i]=false;
        keyrootMap[lml(i)]=i;
    }

    std::map<unsigned int,unsigned int>::const_iterator it;
    for (it=keyrootMap.begin(); it!=keyrootMap.end(); it++) {
        keyroot_[it->second]=true;
    }

}

#endif
