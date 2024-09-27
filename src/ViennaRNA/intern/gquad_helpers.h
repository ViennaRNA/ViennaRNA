#ifndef   VRNA_GQUAD_INTERN_H
#define   VRNA_GQUAD_INTERN_H


#ifndef INLINE
# ifdef __GNUC__
#   define INLINE inline
# else
#   define INLINE
# endif
#endif

/**
 *  Use this macro to loop over each G-quadruplex
 *  delimited by a and b within the subsequence [c,d]
 */
#define FOR_EACH_GQUAD(a, b, c, d)  \
        for ((a) = (d) - VRNA_GQUAD_MIN_BOX_SIZE + 1; (a) >= (c); (a)--) \
        for ((b) = (a) + VRNA_GQUAD_MIN_BOX_SIZE - 1; \
             (b) <= MIN2((d), (a) + VRNA_GQUAD_MAX_BOX_SIZE - 1); \
             (b)++)

#define FOR_EACH_GQUAD_INC(a, b, c, d)  \
        for ((a) = (c); (a) <= (d) - VRNA_GQUAD_MIN_BOX_SIZE + 1; (a)++) \
        for ((b) = (a) + VRNA_GQUAD_MIN_BOX_SIZE - 1; \
             (b) <= MIN2((d), (a) + VRNA_GQUAD_MAX_BOX_SIZE - 1); \
             (b)++)

/**
 *  This macro does almost the same as FOR_EACH_GQUAD() but keeps
 *  the 5' delimiter fixed. 'b' is the 3' delimiter of the gquad,
 *  for gquads within subsequence [a,c] that have 5' delimiter 'a'
 */
#define FOR_EACH_GQUAD_AT(a, b, c)  \
        for ((b) = (a) + VRNA_GQUAD_MIN_BOX_SIZE - 1; \
             (b) <= MIN2((c), (a) + VRNA_GQUAD_MAX_BOX_SIZE - 1); \
             (b)++)


/*
 * Check whether a G-Quadruplex with layer size L and linker lengths l
 * is valid. Perform action r if it doesn't validate
 */
#define CHECK_GQUAD(L, l, r) \
        do { \
          for (size_t i = 0; i < 3; i++) { \
            if ((l)[i] > VRNA_GQUAD_MAX_LINKER_LENGTH) { \
              vrna_log_warning("G-Quadruplex linker length of %u exceeds maximum length of %u", \
                               (unsigned int)(l)[i], \
                               VRNA_GQUAD_MAX_LINKER_LENGTH); \
              r; \
            } \
            if ((l)[i] < VRNA_GQUAD_MIN_LINKER_LENGTH) { \
              vrna_log_warning("G-Quadruplex linker length of %u below minimum length of %u", \
                               (unsigned int)(l)[i], \
                               VRNA_GQUAD_MIN_LINKER_LENGTH); \
              r; \
            } \
          } \
          if ((L) > VRNA_GQUAD_MAX_STACK_SIZE) { \
            vrna_log_warning("G-Quadruplex stack size of %u exceeds maximum stack size of %u", \
                             (L), \
                             VRNA_GQUAD_MAX_STACK_SIZE); \
            r; \
          } \
          if ((L) < VRNA_GQUAD_MIN_STACK_SIZE) { \
            vrna_log_warning("G-Quadruplex stack size of %u below minimum stack size of %u", \
                             (L), \
                             VRNA_GQUAD_MIN_STACK_SIZE); \
            r; \
          } \
        } while (0)


struct gquad_ali_helper {
  const short         **S;
  const unsigned int  **a2s;
  unsigned int        length;
  unsigned int        n_seq;
  vrna_param_t        *P;
  vrna_exp_param_t    *pf;
  unsigned int        L;
  unsigned int        *l;
};



/*
 #################################
 # Useful static functions       #
 #################################
 */

PRIVATE INLINE unsigned int *
get_g_islands_sub(short *S,
                  unsigned int   i,
                  unsigned int   j)
{
  unsigned int x, *gg;

  gg  = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (j - i + 2));
  gg  -= i - 1;

  if (S[j] == 3)
    gg[j] = 1;

  for (x = j - 1; x >= i; x--)
    if (S[x] == 3)
      gg[x] = gg[x + 1] + 1;

  return gg;
}


PRIVATE INLINE unsigned int *
get_g_islands(short *S)
{
  return get_g_islands_sub(S, 1, (unsigned int)S[0]);
}



/**
 *  IMPORTANT:
 *  If you don't know how to use this function, DONT'T USE IT!
 *
 *  The function pointer this function takes as argument is
 *  used for individual calculations with each g-quadruplex
 *  delimited by [i,j].
 *  The function it points to always receives as first 3 arguments
 *  position i, the stack size L and an array l[3] containing the
 *  individual linker sizes.
 *  The remaining 4 (void *) pointers of the callback function receive
 *  the parameters 'data', 'P', 'aux1' and 'aux2' and thus may be
 *  used to pass whatever data you like to.
 *  As the names of those parameters suggest the convention is that
 *  'data' should be used as a pointer where data is stored into,
 *  e.g the MFE or PF and the 'P' parameter should actually be a
 *  'vrna_param_t *' or 'vrna_exp_param_t *' type.
 *  However, what you actually pass obviously depends on the
 *  function the pointer is pointing to.
 *
 *  Although all of this may look like an overkill, it is found
 *  to be almost as fast as implementing g-quadruplex enumeration
 *  in each individual scenario, i.e. code duplication.
 *  Using this function, however, ensures that all g-quadruplex
 *  enumerations are absolutely identical.
 */
/**
 *  We could've also created a macro that loops over all G-quadruplexes
 *  delimited by i and j. However, for the fun of it we use this function
 *  that receives a pointer to a callback function which in turn does the
 *  actual computation for each quadruplex found.
 */
PRIVATE void
process_gquad_enumeration(unsigned int *gg,
                          unsigned int i,
                          unsigned int j,
                          void ( *f )(unsigned int, unsigned int, unsigned int *,
                                      void *, void *, void *, void *),
                          void *data,
                          void *P,
                          void *aux1,
                          void *aux2) VRNA_UNUSED;


PRIVATE void
process_gquad_enumeration(unsigned int *gg,
                          unsigned int i,
                          unsigned int j,
                          void ( *f )(unsigned int, unsigned int, unsigned int *,
                                      void *, void *, void *, void *),
                          void *data,
                          void *P,
                          void *aux1,
                          void *aux2)
{
  unsigned int L, l[3], n, max_linker, maxl0, maxl1;

  n = j - i + 1;

  if ((n >= VRNA_GQUAD_MIN_BOX_SIZE) &&
      (n <= VRNA_GQUAD_MAX_BOX_SIZE)) {
    for (L = MIN2(gg[i], VRNA_GQUAD_MAX_STACK_SIZE);
         L >= VRNA_GQUAD_MIN_STACK_SIZE;
         L--)
      if (gg[j - L + 1] >= L) {
        max_linker = n - 4 * L;
        if ((max_linker >= 3 * VRNA_GQUAD_MIN_LINKER_LENGTH) &&
            (max_linker <= 3 * VRNA_GQUAD_MAX_LINKER_LENGTH)) {
          maxl0 = MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH,
                       max_linker - 2 * VRNA_GQUAD_MIN_LINKER_LENGTH
                       );
          for (l[0] = VRNA_GQUAD_MIN_LINKER_LENGTH;
               l[0] <= maxl0;
               l[0]++)
            if (gg[i + L + l[0]] >= L) {
              maxl1 = MIN2(VRNA_GQUAD_MAX_LINKER_LENGTH,
                           max_linker - l[0] - VRNA_GQUAD_MIN_LINKER_LENGTH
                           );
              for (l[1] = VRNA_GQUAD_MIN_LINKER_LENGTH;
                   l[1] <= maxl1;
                   l[1]++) {
                if (gg[i + 2 * L + l[0] + l[1]] >= L) {
                  l[2] = max_linker - l[0] - l[1];
                  if ((l[2] <= VRNA_GQUAD_MAX_LINKER_LENGTH) &&
                      (gg[i + 3 * L + l[0] + l[1] + l[2]] >= L)) {
                    f(i, L, &(l[0]), data, P, aux1, aux2);
                  }
                }
              }
            }
        }
      }
  }
}


#endif
