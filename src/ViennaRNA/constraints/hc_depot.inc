/* hard constraints depot handling */

struct hc_basepair {
  unsigned int  i;
  unsigned int  j;
  unsigned int  strand_j;
  unsigned char context;
};


struct hc_nuc {
  unsigned int  i;
  int           direction;
  unsigned char context;
  unsigned char nonspec;
};


/* store, for each strand, a list of nucleotide/base pair constraints */
struct vrna_hc_depot_s {
  unsigned int        strands;
  size_t              *up_size;
  size_t              *up_num;
  struct hc_nuc       **up;
  size_t              *bp_size;
  size_t              *bp_num;
  struct hc_basepair  **bp;
};


/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE void
hc_depot_init(vrna_fold_compound_t *fc);


PRIVATE void
hc_depot_store_up(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          strand,
                  unsigned char         context);


PRIVATE void
hc_depot_store_nonspec(vrna_fold_compound_t  *fc,
                       unsigned int          i,
                       unsigned int          strand,
                       int                   d,
                       unsigned char         context);


PRIVATE void
hc_depot_store_bp(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          strand_i,
                  unsigned int          j,
                  unsigned int          strand_j,
                  unsigned int          context);


PRIVATE void
hc_depot_free(vrna_hc_t *hc);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

PRIVATE void
hc_depot_init(vrna_fold_compound_t *fc)
{
  vrna_hc_t *hc = fc->hc;
  if (!hc->depot) {
    hc->depot = (vrna_hc_depot_t *)vrna_alloc(sizeof(vrna_hc_depot_t));

    /*
        by default, we only allocate memory for potential constraints for
        each strand. Missing entries, i.e. NULL pointers for any of the nt
        a constraint may be specified for are considered unconstrained, i.e
        default rules apply
    */
    hc->depot->strands = fc->strands;
    if (fc->strands > 0) {
      hc->depot->up_size = (size_t *)vrna_alloc(sizeof(size_t) * fc->strands);
      hc->depot->up_num  = (size_t *)vrna_alloc(sizeof(size_t) * fc->strands);
      hc->depot->up = (struct hc_nuc **)vrna_alloc(sizeof(struct hc_nuc *) * fc->strands);
      hc->depot->bp_size = (size_t *)vrna_alloc(sizeof(size_t) * fc->strands);
      hc->depot->bp_num  = (size_t *)vrna_alloc(sizeof(size_t) * fc->strands);
      hc->depot->bp = (struct hc_basepair **)vrna_alloc(sizeof(struct hc_basepair *) * fc->strands);
    } else {
      hc->depot->up = NULL;
      hc->depot->bp = NULL;
    }
  }
}


PRIVATE void
hc_depot_store_up(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          strand,
                  unsigned char         context)
{
  vrna_hc_t *hc = fc->hc;

  hc_depot_init(fc);

  if (hc->depot->up_num[strand] == hc->depot->up_size[strand]) {
    hc->depot->up_size[strand] += 32;
    hc->depot->up[strand] = (struct hc_nuc *)vrna_realloc(hc->depot->up[strand],
                                                          sizeof(struct hc_nuc) * (hc->depot->up_size[strand]));
  }

  hc->depot->up[strand][hc->depot->up_num[strand]].i         = i;
  hc->depot->up[strand][hc->depot->up_num[strand]].context   = context;
  hc->depot->up[strand][hc->depot->up_num[strand]].direction = 0;
  hc->depot->up[strand][hc->depot->up_num[strand]].nonspec   = 0;

  hc->depot->up_num[strand]++;
}


PRIVATE void
hc_depot_store_nonspec(vrna_fold_compound_t  *fc,
                       unsigned int          i,
                       unsigned int          strand,
                       int                   d,
                       unsigned char         context)
{
  vrna_hc_t *hc = fc->hc;

  hc_depot_init(fc);

  if (hc->depot->up_num[strand] == hc->depot->up_size[strand]) {
    hc->depot->up_size[strand] += 32;
    hc->depot->up[strand] = (struct hc_nuc *)vrna_realloc(hc->depot->up[strand],
                                                          sizeof(struct hc_nuc) * (hc->depot->up_size[strand]));
  }

  hc->depot->up[strand][hc->depot->up_num[strand]].i         = i;
  hc->depot->up[strand][hc->depot->up_num[strand]].context   = context;
  hc->depot->up[strand][hc->depot->up_num[strand]].direction = d;
  hc->depot->up[strand][hc->depot->up_num[strand]].nonspec   = 1;

  hc->depot->up_num[strand]++;
}


PRIVATE void
hc_depot_store_bp(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          strand_i,
                  unsigned int          j,
                  unsigned int          strand_j,
                  unsigned int          context)
{
  vrna_hc_t           *hc;

  hc_depot_init(fc);

  hc = fc->hc;

  if (hc->depot->bp_num[strand_i] == hc->depot->bp_size[strand_i]) {
    hc->depot->bp_size[strand_i] += 32;
    hc->depot->bp[strand_i] = (struct hc_basepair *)vrna_realloc(hc->depot->bp[strand_i],
                                                               sizeof(struct hc_basepair) * (hc->depot->bp_size[strand_i]));
  }

  hc->depot->bp[strand_i][hc->depot->bp_num[strand_i]].i        = i;
  hc->depot->bp[strand_i][hc->depot->bp_num[strand_i]].j        = j;
  hc->depot->bp[strand_i][hc->depot->bp_num[strand_i]].strand_j = strand_j;
  hc->depot->bp[strand_i][hc->depot->bp_num[strand_i]].context  = context;

  hc->depot->bp_num[strand_i]++;
}


PRIVATE void
hc_depot_free(vrna_hc_t *hc)
{
  unsigned int    s, i, n;
  vrna_hc_depot_t *depot = hc->depot;

  if (depot) {
    if (depot->up) {
      for (s = 0; s < depot->strands; s++)
        free(depot->up[s]);

      free(depot->up);
    }

    if (depot->bp) {
      for (s = 0; s < depot->strands; s++)
        free(depot->bp[s]);

      free(depot->bp);
    }

    free(depot);
  }
  
  hc->depot = NULL;
}

