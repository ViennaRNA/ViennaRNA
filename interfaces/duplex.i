/**********************************************/
/* BEGIN interface for duplex prediction      */
/**********************************************/

%ignore duplexT;

/*
 *  Here, we create a wrapper to support native target language lists
 *  with elements of type duplexT, for which we create a new type
 *  duplex_list_t. By doing so, we can let SWIG handle the std::vector
 *  and calling our destructors for each of the elements when required.
 */
%{

typedef struct {
  int i;
  int j;
  double energy;
  std::string structure;
} duplex_list_t;

%}


typedef struct {
  int i;
  int j;
  double energy;
  std::string structure;
} duplex_list_t;


namespace std {
  %template(DuplexVector) std::vector<duplex_list_t>;
};


%rename (duplexfold) my_duplexfold;
%rename (duplex_subopt) my_duplex_subopt;
%rename (aliduplexfold) my_aliduplexfold;
%rename (aliduplex_subopt) my_aliduplex_subopt;

%{

  duplex_list_t
  my_duplexfold(std::string s1,
                std::string s2)
  {
    duplexT r = duplexfold(s1.c_str(), s2.c_str());
    duplex_list_t res{r.i, r.j, r.energy, r.structure};
    free(r.structure);
    return res;
  }

  std::vector<duplex_list_t>
  my_duplex_subopt( std::string s1,
                    std::string s2,
                    int         delta,
                    int         w)
  {
    std::vector<duplex_list_t> ret;
    duplexT *list, *ptr;
    list = duplex_subopt(s1.c_str(), s2.c_str(), delta, w);
 
    for (ptr = list; ptr->structure != NULL; ptr++) {
      duplex_list_t a{ ptr->i, ptr->j, ptr->energy, ptr->structure};
      ret.push_back(a);
      free(ptr->structure);
    }
    free(list);

    return ret;
  }

  duplex_list_t
  my_aliduplexfold(std::vector<std::string> alignment1,
                   std::vector<std::string> alignment2)
  {
    std::vector<const char*> aln_vec1;
    std::transform(alignment1.begin(), alignment1.end(), std::back_inserter(aln_vec1), convert_vecstring2veccharcp);
    aln_vec1.push_back(NULL); /* mark end of sequences */
    std::vector<const char*> aln_vec2;
    std::transform(alignment2.begin(), alignment2.end(), std::back_inserter(aln_vec2), convert_vecstring2veccharcp);
    aln_vec2.push_back(NULL); /* mark end of sequences */

    duplexT r = aliduplexfold((const char **)&aln_vec1[0], (const char **)&aln_vec2[0]);
    duplex_list_t res{r.i, r.j, r.energy, r.structure};
    free(r.structure);
    return res;

  }

  std::vector<duplex_list_t>
  my_aliduplex_subopt(std::vector<std::string> alignment1,
                      std::vector<std::string> alignment2,
                      int                      delta,
                      int                      w)
  {
    std::vector<duplex_list_t> ret;
    duplexT *list, *ptr;
    std::vector<const char*> aln_vec1;
    std::transform(alignment1.begin(), alignment1.end(), std::back_inserter(aln_vec1), convert_vecstring2veccharcp);
    aln_vec1.push_back(NULL); /* mark end of sequences */
    std::vector<const char*> aln_vec2;
    std::transform(alignment2.begin(), alignment2.end(), std::back_inserter(aln_vec2), convert_vecstring2veccharcp);
    aln_vec2.push_back(NULL); /* mark end of sequences */

    list = aliduplex_subopt((const char **)&aln_vec1[0], (const char **)&aln_vec2[0], delta, w);
    for (ptr = list; ptr->structure != NULL; ptr++) {
      duplex_list_t a{ ptr->i, ptr->j, ptr->energy, ptr->structure};
      ret.push_back(a);
      free(ptr->structure);
    }
    free(list);

    return ret;
  }
%}

#ifdef SWIGPYTHON
%feature("autodoc") my_duplexfold;
%feature("kwargs") my_duplexfold;
%feature("autodoc") my_duplex_subopt;
%feature("kwargs") my_duplex_subopt;
%feature("autodoc") my_aliduplexfold;
%feature("kwargs") my_aliduplexfold;
%feature("autodoc") my_aliduplex_subopt;
%feature("kwargs") my_aliduplex_subopt;
#endif

duplex_list_t my_duplexfold(std::string s1, std::string s2);

std::vector<duplex_list_t> my_duplex_subopt(std::string s1, std::string s2, int delta, int w);

duplex_list_t my_aliduplexfold(std::vector<std::string> alignment1,
                         std::vector<std::string> alignment2);

std::vector<duplex_list_t>  my_aliduplex_subopt(std::vector<std::string> alignment1,
                                                std::vector<std::string> alignment2,
                                                int delta,
                                                int w);

%ignore duplexfold;
%ignore duplex_subopt;
%ignore aliduplexfold;
%ignore aliduplex_subopt;

%include <ViennaRNA/duplex.h>
