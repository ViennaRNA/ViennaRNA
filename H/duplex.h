#ifndef __VIENNA_RNA_PACKAGE_DUPLEX_H__
#define __VIENNA_RNA_PACKAGE_DUPLEX_H__

#include "data_structures.h"

/**
 *  \file duplex.h
 *  \brief Duplex folding function declarations...
 */


duplexT duplexfold( const char *s1,
                    const char *s2);

duplexT *duplex_subopt( const char *s1,
                        const char *s2,
                        int delta,
                        int w);

duplexT aliduplexfold(const char *s1[],
                      const char *s2[]);

duplexT *aliduplex_subopt(const char *s1[],
                          const char *s2[],
                          int delta,
                          int w);

#endif
