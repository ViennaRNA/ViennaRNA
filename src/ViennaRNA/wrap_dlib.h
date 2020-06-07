#ifndef VIENNARNA_DLIB_WRAPPER_H
#define VIENNARNA_DLIB_WRAPPER_H

#ifdef __cplusplus
extern "C" {
#endif

double *
vrna_equilibrium_conc(const double        *eq_constants,
                      double              *concentration_strands,
                      const unsigned int  **A,
                      size_t              num_strands,
                      size_t              num_complexes);


#ifdef __cplusplus
}
#endif

#endif
