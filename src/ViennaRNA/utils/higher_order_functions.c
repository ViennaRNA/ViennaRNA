#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"

#if VRNA_WITH_SIMD_AVX512
#include <immintrin.h>
#elif VRNA_WITH_SIMD_SSE41
#include <emmintrin.h>
#include <smmintrin.h>

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


#endif

PUBLIC int
vrna_fun_zip_add_min(const int  *e1,
                     const int  *e2,
                     int        count)
{
  int     i       = 0;
  int     decomp  = INF;

#if VRNA_WITH_SIMD_AVX512
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
  /* end sources_034_578/modular_decomposition_id3.c */
#elif VRNA_WITH_SIMD_SSE41
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
#endif

  /*
   *  second for loop used by both 128 bit SSE and 512 bit AVX code,
   *  or when compiled without SIMD extensions
   */
  for (; i < count; i++) {
    if ((e1[i] != INF) && (e2[i] != INF)) {
      const int en = e1[i] + e2[i];
      decomp = MIN2(decomp, en);
    }
  }

  return decomp;
}
