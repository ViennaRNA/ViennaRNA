#ifndef PAIR_H_
#define PAIR_H_

#include<list>
#include<iostream>

struct Pair {
  std::list<Pair*> *children;
  unsigned int position;
  Pair(std::list<Pair*> *c, unsigned int p) : children(c), position(p) {}


	void skipit() {
   while (children->size() == 1) {
     std::list<Pair*> t = *children->front()->children;
     *children = t;
   }
   for (std::list<Pair*>::iterator i = children->begin();
        i != children->end(); ++i)
     (*i)->skipit();
	}


  void put(std::list<int> & anchorlist, std::ostream & o) {
    anchorlist.push_back(position);
    for (std::list<Pair*>::iterator i = children->begin(); i != children->end();
         ++i) {
      (*i)->put(anchorlist, o);
    }
  }
};

#endif
