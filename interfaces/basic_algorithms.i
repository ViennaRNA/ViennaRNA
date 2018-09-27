/*#######################################*/
/* Basic algorithms section              */
/*#######################################*/

/**********************************************/
/* BEGIN interface for maximum matching       */
/**********************************************/

%ignore maximumMatching;
%ignore maximumMatchingConstraint;
%ignore maximumMatching2Constraint;

%rename (maximum_matching)  my_maximum_matching;

%{
  int
  my_maximum_matching(std::string sequence)
  {
    return vrna_maximum_matching_simple(sequence.c_str());
  }

%}

int my_maximum_matching(std::string sequence);

%extend vrna_fold_compound_t {

  int
  maxmimum_matching(void)
  {
    return vrna_maximum_matching($self);
  }

}


%include  <ViennaRNA/mm.h>
