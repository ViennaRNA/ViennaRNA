#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"

#include <immintrin.h>


PUBLIC int
vrna_fun_zip_add_min_avx512(const int *e1,
                            const int *e2,
                            int       count)
{
  int     i       = 0;
  int     decomp  = INF;

  __m512i inf = _mm512_set1_epi32(INF);

  /* WBL 21 Aug 2018 Add SSE512 code from sources_034_578/modular_decomposition_id3.c by hand */
  for (i = 0; i < count - 15; i += 16) {
    __m512i   a = _mm512_loadu_si512((__m512i *)&e1[i]);
    __m512i   b = _mm512_loadu_si512((__m512i *)&e2[i]);

    /* compute mask for entries where both, a and b, are less than INF */
    __mmask16 mask = _kand_mask16(_mm512_cmplt_epi32_mask(a, inf),
                                  _mm512_cmplt_epi32_mask(b, inf));

    /* add values */
    __m512i   c = _mm512_add_epi32(a, b);

    /* reduce to minimum (only those where one of the source values was not INF before) */
    const int en = _mm512_mask_reduce_min_epi32(mask, c);

    decomp = MIN2(decomp, en);
  }

  for (; i < count; i++) {
    if ((e1[i] != INF) && (e2[i] != INF)) {
      const int en = e1[i] + e2[i];
      decomp = MIN2(decomp, en);
    }
  }

  return decomp;
}
