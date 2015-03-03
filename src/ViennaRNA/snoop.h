#ifndef VIENNA_RNA_PACKAGE_SNOOP_H
#define VIENNA_RNA_PACKAGE_SNOOP_H

#include <ViennaRNA/data_structures.h>
/** 
*** computes snoRNA-RNA interactions in RNAduplex manner
**/

snoopT snoopfold( const char *s1,
                  const char *s2, 
                  const int penalty,
                  const int threshloop, 
                  const int threshLE,
                  const int threshRE,
                  const int threshDE,
                  const int threshD,
                  const int half_stem,
                  const int max_half_stem,
                  const int min_s2,
                  const int max_s2,
                  const int min_s1,
                  const int max_s1,
                  const int min_d1,
                  const int min_d2,
		  const int fullStemEnergy);

/** 
*** computes snoRNA-RNA suboptimal interactions in RNAduplex manner
**/


snoopT *snoop_subopt( const char *s1,
                      const char *s2,
                      int delta,
                      int w,
                      const int penalty,
                      const int threshloop, 
                      const int threshLE,
                      const int threshRE,
                      const int threshDE,
                      const int threshTE,
                      const int threshSE,
                      const int threshD,
                      const int distance,
                      const int half_stem,
                      const int max_half_stem,
                      const int min_s2,
                      const int max_s2,
                      const int min_s1,
                      const int max_s1,
                      const int min_d1,
                      const int min_d2,
		      const int fullStemEnergy);

/** 
*** computes snoRNA-RNA suboptimal interactions in a RNAplex manner
**/



void Lsnoop_subopt( const char *s1,
                    const char *s2,
                    int delta,
                    int w, 
                    const int penalty,
                    const int threshloop, 
                    const int threshLE,
                    const int threshRE,
                    const int threshDE,
                    const int threshTE,
                    const int threshSE,
                    const int threshD,
                    const int distance,
                    const int half_stem,
                    const int max_half_stem,
                    const int min_s2,
                    const int max_s2,
                    const int min_s1,
                    const int max_s1,
                    const int min_d1,
                    const int min_d2,
                    const int alignment_length,
                    const char* name,
		    const int fullStemEnergy);

/** 
*** computes snoRNA-RNA suboptimal interactions in a RNAplex manner. The stem energy is saved into a list of struct, leading to a runtime improvement of 20%
**/



void Lsnoop_subopt_list ( const char *s1,
                          const char *s2,
                          int delta,
                          int w, 
                          const int penalty,
                          const int threshloop, 
                          const int threshLE,
                          const int threshRE,
                          const int threshDE,
                          const int threshTE,
                          const int threshSE,
                          const int threshD,
                          const int distance,
                          const int half_stem,
                          const int max_half_stem,
                          const int min_s2,
                          const int max_s2,
                          const int min_s1,
                          const int max_s1,
                          const int min_d1,
                          const int min_d2,
                          const int alignment_length,
                          const char *name,
			  const int fullStemEnergy);

/** 
*** computes snoRNA-RNA suboptimal interactions in a RNAplex manner. The stem energy is saved into a list of struct, leading to a runtime improvement of 20%. It considers accessibility
**/


void Lsnoop_subopt_list_XS (const char *s1,
                            const char *s2,
                            const int **access_s1,
                            int delta,
                            int w, 
                            const int penalty,
                            const int threshloop, 
                            const int threshLE,
                            const int threshRE,
                            const int threshDE,
                            const int threshTE,
                            const int threshSE,
                            const int threshD,
                            const int distance,
                            const int half_stem,
                            const int max_half_stem,
                            const int min_s2,
                            const int max_s2,
                            const int min_s1,
                            const int max_s1,
                            const int min_d1,
                            const int min_d2,
                            const int alignment_length,
                            const char *name,
			    const int fullStemEnergy);


/** 
*** computes snoRNA-RNA suboptimal interactions in a RNAduplex manner, and considers accessibility
**/


void snoop_subopt_XS (const char *s1,
                      const char *s2,
                      const int **access_s1,
                      int delta,
                      int w, 
                      const int penalty,
                      const int threshloop, 
                      const int threshLE,
                      const int threshRE,
                      const int threshDE,
                      const int threshTE,
                      const int threshSE,
                      const int threshD,
                      const int distance,
                      const int half_stem,
                      const int max_half_stem,
                      const int min_s2,
                      const int max_s2,
                      const int min_s1,
                      const int max_s1,
                      const int min_d1,
                      const int min_d2,
                      const int alignment_length,
                      const char *name,
		      const int fullStemEnergy);

/**
*** aliduplex-like alignment version of snoop_subopt
 **/

snoopT *alisnoop_subopt(const char **s1,
                        const char **s2,
                        int delta,
                        int w,
                        const int penalty,
                        const int threshloop, 
                        const int threshLE,
                        const int threshRE,
                        const int threshDE,
                        const int threshTE,
                        const int threshSE,
                        const int threshD,
                        const int distance,
                        const int half_stem,
                        const int max_half_stem,
                        const int min_s2,
                        const int max_s2,
                        const int min_s1,
                        const int max_s1,
                        const int min_d1,
                        const int min_d2);

/**
*** RNAplex-like Alignment version of snoop_subopt
 **/



snoopT *aliLsnoop_subopt_list ( const char **s1,
                                const char **s2,
                                int delta,
                                int w, 
                                const int penalty,
                                const int threshloop, 
                                const int threshLE,
                                const int threshRE,
                                const int threshDE,
                                const int threshTE,
                                const int threshSE,
                                const int threshD,
                                const int distance,
                                const int half_stem,
                                const int max_half_stem,
                                const int min_s2,
                                const int max_s2,
                                const int min_s1,
                                const int max_s1,
                                const int min_d1,
                                const int min_d2,
                                const int alignment_length);
/**
*** RNAaliduplex-like version of snoopfold
**/


snoopT alisnoopfold(const char **s1,
                    const char **s2, 
                    const int penalty,
                    const int threshloop,
                    const int threshLE,
                    const int threshRE,
                    const int threshDE,
                    const int threshD,
                    const int half_stem,
                    const int max_half_stem,
                    const int min_s2,
                    const int max_s2,
                    const int min_s1,
                    const int max_s1,
                    const int min_d1,
                    const int min_d2);
/**
*** RNAduplex-like version of snoopfold with accessibility information
**/ 

snoopT snoopfold_XS(const char *s1,
                    const char *s2,
                    const int **access_s1,
                    const int pos,
                    const int max_pos_j,
                    const int penalty,
                    const int threshloop, 
                    const int threshLE,
                    const int threshRE,
                    const int threshDE,
                    const int threshD,
                    const int half_stem,
                    const int max_half_stem,
                    const int min_s2,
                    const int max_s2,
                    const int min_s1,
                    const int max_s1,
                    const int min_d1,
                    const int min_d2,
		    const int fullStemEnergy);




extern int snoop_subopt_sorted;
#endif
