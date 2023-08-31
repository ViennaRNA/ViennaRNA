/*
 * A collection of useful functions to detect various CPU features
 *
 * (c) 2018 - Ronny Lorenz - ViennaRNA Package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/cpu.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/*
 *  cpuid register flags for (64bit) x86 CPUs
 *  See also https://en.wikipedia.org/wiki/CPUID
 */

#define bit_SSE2      (1 << 26) /* stored in EDX after cpuid with EAX=1 */
#define bit_SSE3      (1 << 0)  /* stored in ECX after cpuid with EAX=1 */
#define bit_SSE41     (1 << 19) /* stored in ECX after cpuid with EAX=1 */
#define bit_SSE42     (1 << 20) /* stored in ECX after cpuid with EAX=1 */
#define bit_AVX       (1 << 28) /* stored in ECX after cpuid with EAX=1 */
#define bit_AVX2      (1 << 5)  /* stored in EBX after cpuid with EAX=7, ECX=0 */
#define bit_AVX512F   (1 << 16) /* stored in EBX after cpuid with EAX=7, ECX=0 */

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE INLINE int
execute_cpuid(uint32_t *regs);


PRIVATE unsigned int
cpu_feature_bits(void);


PRIVATE unsigned int
cpu_extended_feature_bits(void);


/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC char *
vrna_cpu_vendor_string(void)
{
  static char name[13] = {
    0
  };
  uint32_t    regs[4] = {
    0
  };

  if (execute_cpuid(&regs[0])) {
    memcpy(name + 0, &regs[1], 4);
    memcpy(name + 4, &regs[3], 4);
    memcpy(name + 8, &regs[2], 4);
    name[12] = '\0';
  }

  return name;
}


PUBLIC unsigned int
vrna_cpu_simd_capabilities(void)
{
  unsigned int capabilities = VRNA_CPU_SIMD_NONE;

  capabilities  |= cpu_feature_bits();
  capabilities  |= cpu_extended_feature_bits();

  return capabilities;
}


/*
 #################################
 # STATIC helper functions below #
 #################################
 */

/*
 * execute 'cpuid' command with values stored in registers regs[0]..regs[3]
 * that directly correspond to EAX, EBX, ECX, and EDX. The result of the
 * 'cpuid' command is then returned in the same register array
 */
PRIVATE INLINE int
execute_cpuid(uint32_t *regs)
{
#if defined(__x86_64__) || defined(_M_AMD64) || defined (_M_X64)
# ifdef _WIN32
#   ifndef __MINGW32__
  __cpuidex(regs, regs[0], regs[2]);
#   else
  __asm__ __volatile__ ("cpuid"
                        : "=a" (regs[0]),
                        "=b" (regs[1]),
                        "=c" (regs[2]),
                        "=d" (regs[3])
                        : "0" (regs[0]),
                        "2" (regs[2]));
#   endif
# else
  __asm__ __volatile__ ("cpuid"
                        : "=a" (regs[0]),
                        "=b" (regs[1]),
                        "=c" (regs[2]),
                        "=d" (regs[3])
                        : "0" (regs[0]),
                        "2" (regs[2]));
# endif
  return 1;
#else
  return 0;
#endif
}


PRIVATE unsigned int
cpu_feature_bits(void)
{
  unsigned int  features = VRNA_CPU_SIMD_NONE;

  uint32_t      regs[4] = {
    1, 0, 0, 0
  };

  if (execute_cpuid(&regs[0])) {
    if (regs[3] & bit_SSE2)
      features |= VRNA_CPU_SIMD_SSE2;

    if (regs[2] & bit_SSE3)
      features |= VRNA_CPU_SIMD_SSE3;

    if (regs[2] & bit_SSE41)
      features |= VRNA_CPU_SIMD_SSE41;

    if (regs[2] & bit_SSE42)
      features |= VRNA_CPU_SIMD_SSE42;

    if (regs[2] & bit_AVX)
      features |= VRNA_CPU_SIMD_AVX;
  }

  return features;
}


PRIVATE unsigned int
cpu_extended_feature_bits(void)
{
  unsigned int  features = VRNA_CPU_SIMD_NONE;

  uint32_t      regs[4] = {
    7, 0, 0, 0
  };

  if (execute_cpuid(&regs[0])) {
    if (regs[1] & bit_AVX2)
      features |= VRNA_CPU_SIMD_AVX2;

    if (regs[1] & bit_AVX512F)
      features |= VRNA_CPU_SIMD_AVX512F;
  }

  return features;
}
