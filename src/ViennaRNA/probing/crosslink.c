
#include "ViennaRNA/datastructures/array.h"
#include "ViennaRNA/math/functions.h"
#include "ViennaRNA/io/utils.h"
#include "ViennaRNA/constraints/soft_iip.h"
#include "ViennaRNA/params/constants.h"
#include "ViennaRNA/probing/crosslink.h"


typedef struct {
  unsigned int  a;          /* start of 5' interval */
  unsigned int  b;          /* end of 5' interval */
  unsigned int  c;          /* start of 3' interval */
  unsigned int  d;          /* end of 3' interval */
  double        p;          /* probability that there is at least one base pair between the intervals */
} interval_data;


PRIVATE vrna_array(interval_data)
parse_crosslink_data(FILE *fp);


PUBLIC int
vrna_probing_data_xlink_file(vrna_fold_compound_t  *fc,
                             const char            *filename,      
                             double                temperature,    
                             double                scaling,        
                             double                limit,          
                             unsigned int          options)        
{
  int     ret = 0;
  double  gamma, max_e, max_p;
  FILE    *fp;

  if ((fc) &&
      (filename)) {
    if (!(fp = fopen(filename, "r"))) {
      vrna_log_warning("Can't open crosslink-data file \"%s\"!",
                       filename);
      return 0;
    }

    double full_scale = scaling * GASCONST * (K0 + temperature) / 1000.;

    if (options & VRNA_PROBING_METHOD_CROSSLINK_LIMIT_PROBABILITY) {
      /* 'limit' argument is maximum probability, so recompute max energy */
      max_p = ((limit >= 0.) && (limit <= 1.)) ? limit : 1.;
      max_e = vrna_math_fun_dbl_logit(max_p,
                                      VRNA_MATH_FUN_LOGIT_BOUNDS_DEFAULT,
                                      VRNA_MATH_FUN_LOGIT_OOB_VALUES_DEFAULT,
                                      VRNA_MATH_FUN_LOGIT_OPTION_DEFAULT);
      if (!isinf(max_e))
        max_e /= full_scale;
    } else {
      /* assume the 'limit' argument is minimum energy for any interval pair */
      max_e = (limit < 0) ? (limit / -full_scale) : 0.;
      max_p = 1. / (1. + exp(-max_e));
    }

    double bounds[2] = {0., max_p};
    double oob_values[2] = {0., max_e};

    vrna_log_debug("iip from file with probability bounds [%g:%g] and corresponding energy bounds [%g:%g] = [%g:%g]",
                   bounds[0], bounds[1],
                   oob_values[0], oob_values[1],
                   -full_scale * oob_values[0], -full_scale * oob_values[1]);

    vrna_array(interval_data) ips = parse_crosslink_data(fp);

    fclose(fp);

    if (ips) {
      ret = 1;

      for (size_t i = 0; i < vrna_array_size(ips); i++) {
        if (options & VRNA_PROBING_METHOD_CROSSLINK_INPUT_ENERGY) {
          gamma = ips[i].p;
          if (gamma < -full_scale * max_e)
            gamma = -full_scale * max_e;
        } else {
          gamma = -full_scale *
                  vrna_math_fun_dbl_logit(ips[i].p,
                                          bounds,
                                          oob_values,
                                          VRNA_MATH_FUN_LOGIT_OPTION_DEFAULT);
        }

        vrna_log_debug("adding xlink interval [%u:%u] [%u:%u] e = %g",
                       ips[i].a,
                       ips[i].b,
                       ips[i].c,
                       ips[i].d,
                       gamma);

        ret &= vrna_sc_add_iip(fc,
                               ips[i].a,
                               ips[i].b,
                               ips[i].c,
                               ips[i].d,
                               gamma,
                               0);
      }
    }

    vrna_array_free(ips);
  }

  return ret;
}


PRIVATE vrna_array(interval_data)
parse_crosslink_data(FILE *fp)
{
  char    *line;
  int     a, b, c, d;
  double  p;

  vrna_array(interval_data) ips;

  vrna_array_init_size(ips, 42);

  while ((line = vrna_read_line(fp))) {
    if (sscanf(line, "%d %d %d %d %lf", &a, &b, &c, &d, &p) == 5) {
      interval_data r = {
        .a  = a,
        .b  = b,
        .c  = c,
        .d  = d,
        .p  = p
      };

      vrna_array_append(ips, r);
    }

    free(line);
  }

  if (vrna_array_size(ips) == 0) {
    vrna_array_free(ips);
    ips = NULL;
  }

  return ips;
}


