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

#define bit_SSE41     (1 << 19) /* stored in ECX after cpuid with EAX=1 */
#define bit_AVX512F   (1 << 16) /* stored in EBX after cpuid with EAX=7, ECX=0 */

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE INLINE int
execute_cpuid(uint32_t *regs);


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


PUBLIC int
vrna_cpu_sse41(void)
{
  uint32_t regs[4] = {
    1, 0, 0, 0
  };

  if ((execute_cpuid(&regs[0])) &&
      (regs[2] & bit_SSE41))
    return 1;

  return 0;
}


PUBLIC int
vrna_cpu_avx512f(void)
{
  uint32_t regs[4] = {
    7, 0, 0, 0
  };

  if ((execute_cpuid(&regs[0])) &&
      (regs[1] & bit_AVX512F))
    return 1;

  return 0;
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
  __asm__ __volatile__ ("cpuid"
                        : "=a" (regs[0]),
                        "=b" (regs[1]),
                        "=c" (regs[2]),
                        "=d" (regs[3])
                        : "0" (regs[0]),
                        "2" (regs[2]));

  return 1;
#else
  return 0;
#endif
}
