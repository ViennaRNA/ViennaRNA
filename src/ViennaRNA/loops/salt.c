#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <math.h>
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/loops/salt.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif


#define Rods_dist 20.
#define PI 3.141592653589793

extern double kn(int, double);


PRIVATE INLINE double
bjerrum_length(double T);


PRIVATE INLINE double
ionic_strength(double rho);


PRIVATE INLINE double
kappa_inv(double rho, double T);


PRIVATE INLINE double
kappa(double rho, double T);


PRIVATE INLINE double
tau_ds(double T);


PRIVATE INLINE double
tau_ss(double T);


PRIVATE INLINE double
pairing_salt_const(double rho, double T);


PRIVATE double
exp_one(double x);


PRIVATE double
approx_hyper(double y)
{
	double a, b, c;

	a = 1/(pow(y,6.)/pow(2*PI,6.)+1);
	b = pow(y, 4.)/(36*pow(PI, 4.)) - pow(y, 3.)/(24*PI*PI) + y*y/(2*PI*PI) - y/2;
	c = log(2*PI/y) - 1.96351;

	return a*b + (1-a)*c;
};


PRIVATE double
loop_salt_aux(double kmlss, int m, double T)
{
	double a, b;

	a =  Gaz_const_kcal * T * bjerrum_length(T) * m * Backbone_length * tau_ss(T) * tau_ss(T);
	b = log(kmlss) - log(PI/2) + Eular_const + approx_hyper(kmlss) + 1/kmlss * (1 - exp(-kmlss) + kmlss*exp_one(kmlss));

	return a*b*100;
};

PUBLIC int
vrna_salt_loop(int m, double rho, double T)
{
	double correction, kmlss, kmlss_ref;
	
	kmlss_ref = kappa(1., T) * m * Backbone_length;
	kmlss = kappa(rho, T) * m * Backbone_length;

	correction = loop_salt_aux(kmlss, m, T) - loop_salt_aux(kmlss_ref, m, T);

	return (int) (correction + 0.5 - (correction<0)); /* round salt correction*/
};

PUBLIC int
vrna_salt_stack(double rho, double T)
{
	double correction, kn_ref, kappa_ref;

	kappa_ref = kappa(1., T);
	kn_ref = kn(0, Rods_dist*kappa_ref);
	correction = 100*pairing_salt_const(rho, T)*(kn(0, Rods_dist*kappa(rho, T)) - kn_ref);
	printf("correction %lf\n", correction);
	return (int) (correction + 0.5 - (correction<0)); /* round salt correction*/
}
