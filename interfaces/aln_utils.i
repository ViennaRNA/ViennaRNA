/*#######################################*/
/* Interface for alignment utilities     */
/*#######################################*/

%rename (aln_consensus_sequence) my_consensus_sequence;
%rename (consensus) my_consensus_sequence;
%rename (aln_consensus_mis) my_aln_consensus_mis;
%rename (consens_mis) my_aln_consensus_mis;

%{
#include <vector>

  std::string
  my_consensus_sequence(std::vector<std::string>  alignment,
                        vrna_md_t                 *md_p = NULL)
  {
    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  v;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(v), convert_vecstring2veccharcp);
    v.push_back(NULL); /* mark end of sequences */

    char *c = vrna_aln_consensus_sequence((const char **)&v[0], md_p);
    std::string cons(c);
    free(c);
    return cons;
  }

  std::string
  my_aln_consensus_mis(std::vector<std::string> alignment,
                       vrna_md_t                *md_p = NULL)
  {
    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  v;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(v), convert_vecstring2veccharcp);
    v.push_back(NULL); /* mark end of sequences */

    char *c = vrna_aln_consensus_mis((const char **)&v[0], md_p);
    std::string mis(c);
    free(c);
    return mis;
  }

%}

#ifdef SWIGPYTHON
%feature("autodoc") my_consensus_sequence;
%feature("kwargs") my_consensus_sequence;
%feature("autodoc") my_aln_consensus_mis;
%feature("kwargs") my_aln_consensus_mis;
#endif

std::string
my_consensus_sequence(std::vector<std::string> alignment, vrna_md_t *md_p = NULL);

std::string
my_aln_consensus_mis(std::vector<std::string> alignment, vrna_md_t *md_p = NULL);

%ignore consensus;
%ignore consens_mis;


%rename (aln_mpi) my_aln_mpi;
%{
#include <vector>

  int
  my_aln_mpi(std::vector<std::string> alignment)
  {
    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  v;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(v), convert_vecstring2veccharcp);
    v.push_back(NULL); /* mark end of sequences */

    int mpi = vrna_aln_mpi((const char **)&v[0]);

    return mpi;
  }

%}

#ifdef SWIGPYTHON
%feature("autodoc") my_aln_mpi;
%feature("kwargs") my_aln_mpi;
#endif

int
my_aln_mpi(std::vector<std::string> alignment);


%rename (aln_pscore) my_aln_pscore;

%{
#include <vector>

  std::vector<std::vector<int> >
  my_aln_pscore(std::vector<std::string>  alignment,
                vrna_md_t                 *md = NULL)
  {

    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  v;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(v), convert_vecstring2veccharcp);
    v.push_back(NULL); /* mark end of sequences */

    std::vector<std::vector<int> > pscore;
    int *ps = vrna_aln_pscore((const char **)&v[0], md);

    int n     = alignment[0].length();
    int *idx  = vrna_idx_col_wise(n);

    std::vector<int> z_row(n+1, 0);
    pscore.push_back(z_row);

    for(int i = 1; i < n; i++){
      std::vector<int> score_i;
      score_i.push_back(0);
      for(int j = 1; j <= i; j++)
        score_i.push_back(ps[idx[i] + j]);
      for(int j = i + 1; j <= n; j++)
        score_i.push_back(ps[idx[j] + i]);
      pscore.push_back(score_i);
    }

    free(ps);
    free(idx);

    return pscore;
  }

%}

#ifdef SWIGPYTHON
%feature("autodoc") my_aln_pscore;
%feature("kwargs") my_aln_pscore;
#endif

std::vector<std::vector<int> >
my_aln_pscore(std::vector<std::string> alignment,
              vrna_md_t *md = NULL);


%rename (aln_conservation_struct) my_aln_conservation_struct;
%rename (aln_conservation_col) my_aln_conservation_col;

%{
#include <vector>

  std::vector<double>
  my_aln_conservation_struct(std::vector<std::string> alignment,
                             std::string              structure,
                             vrna_md_t                *md = NULL)
  {
    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  v;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(v), convert_vecstring2veccharcp);
    v.push_back(NULL); /* mark end of sequences */

    std::vector<double> conservation;

    float *c = vrna_aln_conservation_struct((const char **)&v[0], structure.c_str(), md);

    if (c) {
      for (unsigned int i = 0; i <= alignment[0].size(); i++)
        conservation.push_back((double)c[i]);

      free(c);
    }

    return conservation;
  }
  std::vector<double>
  my_aln_conservation_col(std::vector<std::string> alignment,
                          vrna_md_t                *md = NULL,
                          unsigned int             options = VRNA_MEASURE_SHANNON_ENTROPY)
  {
    /* convert std::vector<std::string> to vector<const char *> */
    std::vector<const char*>  v;
    std::transform(alignment.begin(), alignment.end(), std::back_inserter(v), convert_vecstring2veccharcp);
    v.push_back(NULL); /* mark end of sequences */

    std::vector<double> conservation;

    float *c = vrna_aln_conservation_col((const char **)&v[0], md, options);

    if (c) {
      for (unsigned int i = 0; i <= alignment[0].size(); i++)
        conservation.push_back((double)c[i]);

      free(c);
    }

    return conservation;
  }
%}


#ifdef SWIGPYTHON
%feature("autodoc") my_aln_conservation_struct;
%feature("autodoc") my_aln_conservation_col;
%feature("kwargs") my_aln_conservation_struct;
%feature("kwargs") my_aln_conservation_col;
#endif

std::vector<double>
my_aln_conservation_struct(std::vector<std::string> alignment,
                           std::string structure,
                           vrna_md_t *md = NULL);

std::vector<double>
my_aln_conservation_col(std::vector<std::string> alignment,
                        vrna_md_t                *md = NULL,
                        unsigned int             options = VRNA_MEASURE_SHANNON_ENTROPY);


%ignore read_clustal;
%ignore get_ungapped_sequence;
%ignore get_mpi;
%ignore encode_ali_sequence;
%ignore alloc_sequence_arrays;
%ignore free_sequence_arrays;

%constant unsigned int ALN_DEFAULT              = VRNA_ALN_DEFAULT;
%constant unsigned int ALN_RNA                  = VRNA_ALN_RNA;
%constant unsigned int ALN_DNA                  = VRNA_ALN_DNA;
%constant unsigned int ALN_UPPERCASE            = VRNA_ALN_UPPERCASE;
%constant unsigned int ALN_LOWERCASE            = VRNA_ALN_LOWERCASE;
%constant unsigned int MEASURE_SHANNON_ENTROPY  = VRNA_MEASURE_SHANNON_ENTROPY;

%include  <ViennaRNA/utils/alignments.h>
