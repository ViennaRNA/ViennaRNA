#ifndef VIENNA_RNA_PACKAGE_UTILS_CPU_H
#define VIENNA_RNA_PACKAGE_UTILS_CPU_H

#define VRNA_CPU_SIMD_NONE      0U
#define VRNA_CPU_SIMD_SSE2      1U
#define VRNA_CPU_SIMD_SSE3      2U
#define VRNA_CPU_SIMD_SSE41     4U
#define VRNA_CPU_SIMD_SSE42     8U
#define VRNA_CPU_SIMD_AVX       16U
#define VRNA_CPU_SIMD_AVX2      32U
#define VRNA_CPU_SIMD_AVX512F   64U


char *
vrna_cpu_vendor_string(void);


unsigned int
vrna_cpu_simd_capabilities(void);


#endif
