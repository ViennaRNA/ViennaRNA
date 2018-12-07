#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"

#include <emmintrin.h>
#include <smmintrin.h>

static int
horizontal_min_Vec4i(__m128i x);


PUBLIC int
vrna_fun_zip_add_min_sse41(const int  *e1,
                           const int  *e2,
                           int        count)
{
  int     i       = 0;
  int     decomp  = INF;

  __m128i inf = _mm_set1_epi32(INF);

  for (i = 0; i < count - 3; i += 4) {
    __m128i a = _mm_loadu_si128((__m128i *)&e1[i]);
    __m128i b = _mm_loadu_si128((__m128i *)&e2[i]);
    __m128i c = _mm_add_epi32(a, b);

    /* create mask for non-INF values */
    __m128i mask = _mm_and_si128(_mm_cmplt_epi32(a, inf),
                                 _mm_cmplt_epi32(b, inf));

    /* delete results where a or b has been INF before */
    c = _mm_and_si128(mask, c);

    /* fill all values with INF if they've been INF in a or b before */
    __m128i   res = _mm_or_si128(c, _mm_andnot_si128(mask, inf));
    const int en  = horizontal_min_Vec4i(res);

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


/*
 *  SSE minimum
 *  see also: http://stackoverflow.com/questions/9877700/getting-max-value-in-a-m128i-vector-with-sse
 */
static int
horizontal_min_Vec4i(__m128i x)
{
  __m128i min1  = _mm_shuffle_epi32(x, _MM_SHUFFLE(0, 0, 3, 2));
  __m128i min2  = _mm_min_epi32(x, min1);
  __m128i min3  = _mm_shuffle_epi32(min2, _MM_SHUFFLE(0, 0, 0, 1));
  __m128i min4  = _mm_min_epi32(min2, min3);

  return _mm_cvtsi128_si32(min4);
}
