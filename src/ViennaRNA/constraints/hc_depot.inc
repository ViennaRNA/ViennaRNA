/* hard constraints depot handling */

struct hc_basepair {
  size_t        list_size;
  size_t        list_mem;
  unsigned int  *j;
  unsigned int  *strand_j;
  unsigned char *context;
};

struct hc_nuc {
  int           direction;
  unsigned char context;
  unsigned char nonspec;
};


/* store, for each strand, a list of nucleotide/base pair constraints */
struct vrna_hc_depot_s {
  unsigned int            strands;
  size_t                  *up_size;
  struct hc_nuc           **up;
  size_t                  *bp_size;
  struct hc_basepair      **bp;
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
                  unsigned char          context);


PRIVATE void
hc_depot_free(vrna_hc_t *hc);


PRIVATE void
hc_depot_resize_bp(struct vrna_hc_depot_s *depot,
                   unsigned int           strand,
                   unsigned int           i);


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
      hc->depot->up = (struct hc_nuc **)vrna_alloc(sizeof(struct hc_nuc *) * fc->strands);
      hc->depot->bp_size = (size_t *)vrna_alloc(sizeof(size_t) * fc->strands);
      hc->depot->bp = (struct hc_basepair **)vrna_alloc(sizeof(struct hc_basepair *) * fc->strands);
    } else {
      hc->depot->up_size = NULL;
      hc->depot->up = NULL;
      hc->depot->bp_size = NULL;
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
  size_t    k, old_size;
  vrna_hc_t *hc = fc->hc;

  hc_depot_init(fc);

  if (hc->depot->up_size[strand] < i) {
    old_size = hc->depot->up_size[strand];
    hc->depot->up_size[strand] = i;
    hc->depot->up[strand] = (struct hc_nuc *)vrna_realloc(hc->depot->up[strand],
                                                          sizeof(struct hc_nuc) * (hc->depot->up_size[strand] + 1));
    /* initialize new entries */
    for (k = old_size + 1; k < i; k++) {
      hc->depot->up[strand][k].context   = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS | VRNA_CONSTRAINT_CONTEXT_NO_REMOVE;
      hc->depot->up[strand][k].direction = 0;
      hc->depot->up[strand][k].nonspec   = 0;
    }
  }

  hc->depot->up[strand][i].context   = context;
  hc->depot->up[strand][i].direction = 0;
  hc->depot->up[strand][i].nonspec   = 0;
}


PRIVATE void
hc_depot_store_nonspec(vrna_fold_compound_t  *fc,
                       unsigned int          i,
                       unsigned int          strand,
                       int                   d,
                       unsigned char         context)
{
  size_t    k, old_size;
  vrna_hc_t *hc = fc->hc;

  hc_depot_init(fc);

  if (hc->depot->up_size[strand] < i) {
    old_size = hc->depot->up_size[strand];
    hc->depot->up_size[strand] = i;
    hc->depot->up[strand] = (struct hc_nuc *)vrna_realloc(hc->depot->up[strand],
                                                          sizeof(struct hc_nuc) *
                                                          (hc->depot->up_size[strand] + 1));
    /* initialize new entries */
    for (k = old_size + 1; k < i; k++) {
      hc->depot->up[strand][k].context   = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS | VRNA_CONSTRAINT_CONTEXT_NO_REMOVE;
      hc->depot->up[strand][k].direction = 0;
      hc->depot->up[strand][k].nonspec   = 0;
    }
  }

  hc->depot->up[strand][i].context   = context;
  hc->depot->up[strand][i].direction = d;
  hc->depot->up[strand][i].nonspec   = 1;
}


PRIVATE void
hc_depot_resize_bp(struct vrna_hc_depot_s *depot,
                   unsigned int           strand,
                   unsigned int           i)
{
  size_t  old_size, k;

  /* 1a, resize memory for constraints on current strand if necessary */
  if (depot->bp_size[strand] < i) {
    old_size = depot->bp_size[strand];
    depot->bp_size[strand] = i;
    depot->bp[strand] = (struct hc_basepair *)vrna_realloc(depot->bp[strand],
                                                                 sizeof(struct hc_basepair) *
                                                                 (depot->bp_size[strand] + 1));

    /* initialize new entries */
    for (k = old_size + 1; k <= i; k++) {
      depot->bp[strand][k].list_mem   = 0;
      depot->bp[strand][k].list_size  = 0;
      depot->bp[strand][k].j          = NULL;
      depot->bp[strand][k].strand_j   = NULL;
      depot->bp[strand][k].context    = NULL;
    }
  }

  /* 1b, resize memory for constraints for nucleotide i if necessary */
  if (depot->bp[strand][i].list_size == depot->bp[strand][i].list_mem) {
    depot->bp[strand][i].list_mem += 32;
    depot->bp[strand][i].j = (unsigned int *)vrna_realloc(depot->bp[strand][i].j,
                                                          sizeof(unsigned int) *
                                                          (depot->bp[strand][i].list_mem + 1));
    depot->bp[strand][i].strand_j = (unsigned int *)vrna_realloc(depot->bp[strand][i].strand_j,
                                                                 sizeof(unsigned int) *
                                                                 (depot->bp[strand][i].list_mem + 1));
    depot->bp[strand][i].context = (unsigned char *)vrna_realloc(depot->bp[strand][i].context,
                                                                sizeof(unsigned char) *
                                                                (depot->bp[strand][i].list_mem + 1));
  }
}


PRIVATE void
hc_depot_store_bp(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          strand_i,
                  unsigned int          j,
                  unsigned int          strand_j,
                  unsigned char         context)
{
  size_t    next_entry;
  vrna_hc_t *hc;

  hc_depot_init(fc);

  hc = fc->hc;

  /* 1st, store the actual constraint */
  hc_depot_resize_bp(hc->depot, strand_i, i);

  next_entry = hc->depot->bp[strand_i][i].list_size;

  hc->depot->bp[strand_i][i].j[next_entry]         = j;
  hc->depot->bp[strand_i][i].strand_j[next_entry]  = strand_j;
  hc->depot->bp[strand_i][i].context[next_entry]   = context;

  hc->depot->bp[strand_i][i].list_size++;

  /* 2nd. store the 'reverse-lookup' for this constraint */
  hc_depot_resize_bp(hc->depot, strand_j, j);

  next_entry = hc->depot->bp[strand_j][j].list_size;

  hc->depot->bp[strand_j][j].j[next_entry]         = i;
  hc->depot->bp[strand_j][j].strand_j[next_entry]  = strand_i;
  hc->depot->bp[strand_j][j].context[next_entry]   = context;

  hc->depot->bp[strand_j][j].list_size++;
}


PRIVATE void
hc_depot_free(vrna_hc_t *hc)
{
  unsigned int    s, i;
  vrna_hc_depot_t *depot = hc->depot;

  if (depot) {
    if (depot->up) {
      for (s = 0; s < depot->strands; s++)
        free(depot->up[s]);

      free(depot->up);
    }

    if (depot->bp) {
      for (s = 0; s < depot->strands; s++) {
        for (i = 1; i <= depot->bp_size[s]; i++) {
          free(depot->bp[s][i].j);
          free(depot->bp[s][i].strand_j);
          free(depot->bp[s][i].context);
        }
        free(depot->bp[s]);
      }

      free(depot->bp);
    }

    free(depot->bp_size);
    free(depot->up_size);
    free(depot);
  }
  
  hc->depot = NULL;
}

