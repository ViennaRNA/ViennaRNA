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

struct hc_basepair_rev {
  size_t    source_size;
  size_t    source_mem;
  unsigned int  *strand_i;
  unsigned int  *i;
  size_t        *entry;
};


/* store, for each strand, a list of nucleotide/base pair constraints */
struct vrna_hc_depot_s {
  unsigned int            strands;
  size_t                  *up_size;
  struct hc_nuc           **up;
  size_t                  *bp_size;
  struct hc_basepair      **bp;
  size_t                  *rev_size;
  struct hc_basepair_rev  **rev;
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
      hc->depot->rev_size = (size_t *)vrna_alloc(sizeof(size_t) * fc->strands);
      hc->depot->rev = (struct hc_basepair_rev **)vrna_alloc(sizeof(struct hc_basepair_rev *) * fc->strands);
    } else {
      hc->depot->up_size = NULL;
      hc->depot->up = NULL;
      hc->depot->bp_size = NULL;
      hc->depot->bp = NULL;
      hc->depot->rev_size = NULL;
      hc->depot->rev = NULL;
    }
  }
}


PRIVATE void
hc_depot_store_up(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          strand,
                  unsigned char         context)
{
  unsigned int  k, old_size;
  vrna_hc_t     *hc = fc->hc;

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
hc_depot_store_bp(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          strand_i,
                  unsigned int          j,
                  unsigned int          strand_j,
                  unsigned char         context)
{
  size_t    k, old_size, next_entry, entry;
  vrna_hc_t           *hc;

  hc_depot_init(fc);

  hc = fc->hc;

  /* 1st, store the actual constraint */

  /* 1a, resize memory for constraints on strand_i if necessary */
  if (hc->depot->bp_size[strand_i] < i) {
    old_size = hc->depot->bp_size[strand_i];
    hc->depot->bp_size[strand_i] = i;
    hc->depot->bp[strand_i] = (struct hc_basepair *)vrna_realloc(hc->depot->bp[strand_i],
                                                                 sizeof(struct hc_basepair) *
                                                                 (hc->depot->bp_size[strand_i] + 1));
    /* initialize new entries */
    for (k = old_size + 1; k <= i; k++) {
      hc->depot->bp[strand_i][k].list_mem   = 0;
      hc->depot->bp[strand_i][k].list_size  = 0;
      hc->depot->bp[strand_i][k].j          = NULL;
      hc->depot->bp[strand_i][k].strand_j   = NULL;
      hc->depot->bp[strand_i][k].context    = NULL;
    }
  }

  /* 1b, resize memory for constraints for nucleotide i on strand_i if necessary */
  if (hc->depot->bp[strand_i][i].list_size == hc->depot->bp[strand_i][i].list_mem) {
    hc->depot->bp[strand_i][i].list_mem += 32;
    hc->depot->bp[strand_i][i].j = (unsigned int *)vrna_realloc(hc->depot->bp[strand_i][i].j,
                                                                sizeof(unsigned int) *
                                                                (hc->depot->bp[strand_i][i].list_mem + 1));
    hc->depot->bp[strand_i][i].strand_j = (unsigned int *)vrna_realloc(hc->depot->bp[strand_i][i].strand_j,
                                                                       sizeof(unsigned int) *
                                                                       (hc->depot->bp[strand_i][i].list_mem + 1));
    hc->depot->bp[strand_i][i].context = (unsigned char *)vrna_realloc(hc->depot->bp[strand_i][i].context,
                                                                sizeof(unsigned char) *
                                                                (hc->depot->bp[strand_i][i].list_mem + 1));
  }

  /* 1c, store the constraint */
  entry = next_entry = hc->depot->bp[strand_i][i].list_size;

  hc->depot->bp[strand_i][i].j[next_entry]         = j;
  hc->depot->bp[strand_i][i].strand_j[next_entry]  = strand_j;
  hc->depot->bp[strand_i][i].context[next_entry]   = context;

  hc->depot->bp[strand_i][i].list_size++;

  /* 2nd. store a reverse-lookup entry for this base pair constraint */

  /* 2a, resize memory for reverse lookup entries affecting strand_j if necessary */
  if (hc->depot->rev_size[strand_j] < j) {
    old_size = hc->depot->rev_size[strand_j];
    hc->depot->rev_size[strand_j] = j;
    hc->depot->rev[strand_j] = (struct hc_basepair_rev *)vrna_realloc(hc->depot->rev[strand_j],
                                                                      sizeof(struct hc_basepair_rev) *
                                                                      (hc->depot->rev_size[strand_j] + 1));

    /* initialize new entries */
    for (k = old_size + 1; k <= j; k++) {
      hc->depot->rev[strand_j][k].source_size = 0;
      hc->depot->rev[strand_j][k].source_mem  = 0;
      hc->depot->rev[strand_j][k].strand_i    = NULL;
      hc->depot->rev[strand_j][k].i           = NULL;
      hc->depot->rev[strand_j][k].entry       = NULL;
    }
  }

  /* 2b, resize memory for reverse lookup entries affecting position j on strand_j if necessary */
  if (hc->depot->rev[strand_j][j].source_size == hc->depot->rev[strand_j][j].source_mem) {
    hc->depot->rev[strand_j][j].source_mem += 32;
    hc->depot->rev[strand_j][j].strand_i = (unsigned int *)vrna_realloc(hc->depot->rev[strand_j][j].strand_i,
                                                                             sizeof(unsigned int) *
                                                                             (hc->depot->rev[strand_j][j].source_mem + 1));
    hc->depot->rev[strand_j][j].i = (unsigned int *)vrna_realloc(hc->depot->rev[strand_j][j].i,
                                                                sizeof(unsigned int) *
                                                                (hc->depot->rev[strand_j][j].source_mem + 1));
    hc->depot->rev[strand_j][j].entry = (size_t *)vrna_realloc(hc->depot->rev[strand_j][j].entry,
                                                                sizeof(size_t) *
                                                                (hc->depot->rev[strand_j][j].source_mem + 1));
  }

  /* 2c, store the reverse lookup */
  next_entry = hc->depot->rev[strand_j][j].source_size;

  hc->depot->rev[strand_j][j].strand_i[next_entry]  = strand_i;
  hc->depot->rev[strand_j][j].i[next_entry]         = i;
  hc->depot->rev[strand_j][j].entry[next_entry]     = entry;
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

    if (depot->rev) {
      for (s = 0; s < depot->strands; s++) {
        for (i = 1; i <= depot->rev_size[s]; i++) {
          free(depot->rev[s][i].strand_i);
          free(depot->rev[s][i].i);
          free(depot->rev[s][i].entry);
        }

        free(depot->rev[s]);
      }

      free(depot->rev);
    }

    free(depot);
  }
  
  hc->depot = NULL;
}

