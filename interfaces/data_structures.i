/**********************************************/
/* BEGIN interface for data structures        */
/**********************************************/

/* ignore all data structures we handle elsewhere */
%ignore PAIR;
%ignore plist;
%ignore cpair;
%ignore sect;
%ignore bondT;
%ignore pu_contrib;
%ignore interact;
%ignore pu_out;
%ignore constrain;
%ignore folden;
%ignore snoopT;
%ignore dupVar;
%ignore folden;
%ignore node;
%ignore snoopT;
%ignore dupVar;

%rename(basepair) vrna_basepair_t;

typedef struct {
  int i;
  int j;
} vrna_basepair_t;

#ifdef SWIGPERL5
%rename(_next) node::next;
#endif

%include <ViennaRNA/datastructures/basic.h>
