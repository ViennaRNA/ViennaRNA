#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/cpu.h"
#include "ViennaRNA/utils/higher_order_functions.h"


typedef int (*proto_fun_zip_reduce)(const int  *a,
                                    const int  *b,
                                    int        size);


/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
static int
zip_add_min_dispatcher(const int  *a,
                       const int  *b,
                       int        size);


static int
fun_zip_add_min_default(const int *e1,
                        const int *e2,
                        int       count);


#if VRNA_WITH_SIMD_AVX512
extern int
vrna_fun_zip_add_min_avx512(const int *e1,
                            const int *e2,
                            int       count);


#endif

#if VRNA_WITH_SIMD_SSE41
extern int
vrna_fun_zip_add_min_sse41(const int  *e1,
                           const int  *e2,
                           int        count);


#endif


static proto_fun_zip_reduce fun_zip_add_min = &zip_add_min_dispatcher;


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void
vrna_fun_dispatch_disable(void)
{
  fun_zip_add_min = &fun_zip_add_min_default;
}


PUBLIC void
vrna_fun_dispatch_enable(void)
{
  fun_zip_add_min = &zip_add_min_dispatcher;
}


PUBLIC int
vrna_fun_zip_add_min(const int  *e1,
                     const int  *e2,
                     int        count)
{
  return (*fun_zip_add_min)(e1, e2, count);
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */

/* zip_add_min() dispatcher */
static int
zip_add_min_dispatcher(const int  *a,
                       const int  *b,
                       int        size)
{
  unsigned int features = vrna_cpu_simd_capabilities();

#if VRNA_WITH_SIMD_AVX512
  if (features & VRNA_CPU_SIMD_AVX512F) {
    fun_zip_add_min = &vrna_fun_zip_add_min_avx512;
    goto exec_fun_zip_add_min;
  }

#endif

#if VRNA_WITH_SIMD_SSE41
  if (features & VRNA_CPU_SIMD_SSE41) {
    fun_zip_add_min = &vrna_fun_zip_add_min_sse41;
    goto exec_fun_zip_add_min;
  }

#endif

  fun_zip_add_min = &fun_zip_add_min_default;

exec_fun_zip_add_min:

  return (*fun_zip_add_min)(a, b, size);
}


static int
fun_zip_add_min_default(const int *e1,
                        const int *e2,
                        int       count)
{
  int i;
  int decomp = INF;

  for (i = 0; i < count; i++) {
    if ((e1[i] != INF) && (e2[i] != INF)) {
      const int en = e1[i] + e2[i];
      decomp = MIN2(decomp, en);
    }
  }

  return decomp;
}
