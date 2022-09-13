#ifndef VIENNA_RNA_PACKAGE_LOOPS_SALT_H
#define VIENNA_RNA_PACKAGE_LOOPS_SALT_H


#include <math.h>
#include "ViennaRNA/utils/basic.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#define Helical_Rise        3.4
#define Helical_Rise_inv    0.2941
#define Gaz_const_kcal      0.001987
#define Backbone_length     6.4
#define Eular_const         0.58


PRIVATE INLINE double
bjerrum_length(double T)
{
  return 2088.762/T;
};


PRIVATE INLINE double
ionic_strength(double rho)
{
  return rho;
};


PRIVATE INLINE double
kappa_inv(double rho, double T)
{
  return 8.1285/(sqrt(bjerrum_length(T)*ionic_strength(rho)));
};


PRIVATE INLINE double
kappa(double rho, double T)
{
  /* return sqrt(bjerrum_length(T)*ionic_strength(rho))/8.1285; */
  return 1/kappa_inv(rho, T);
};


PRIVATE INLINE double
tau_ds(double T)
{
  double bjerrum_length_inv;

  bjerrum_length_inv = 1/bjerrum_length(T);
  return MIN2(Helical_Rise_inv, bjerrum_length_inv);
};


PRIVATE INLINE double
tau_ss(double T)
{
  double bjerrum_length_inv;

  bjerrum_length_inv = 1/bjerrum_length(T);
  return MIN2(1/Backbone_length, bjerrum_length_inv);
};


PRIVATE INLINE double
pairing_salt_const(double rho, double T)
{
  return 2*Gaz_const_kcal*T*bjerrum_length(T)*Helical_Rise*tau_ds(T)*tau_ds(T);
};


/* TODO: Define incomplete gamma function */
PRIVATE double
exp_one(double x)
{
	return 1.;
};


int
vrna_salt_loop(int m, double salt, double T);

int
vrna_salt_stack(double salt, double T);


#endif
