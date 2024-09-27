#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <math.h>
#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/model.h"
#include "ViennaRNA/params/constants.h"
#include "ViennaRNA/params/salt.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif


#define Rods_dist 20.
#define PI 3.141592653589793
/* #define Backbone_length     6.4 */
#define Eular_const         0.58

#define roundint(x) ((int) (x + 0.5 - (x<0)))


extern double kn(int, double);
extern double expn(int, double);

/* Temperature-dependent epsilonr */
PRIVATE INLINE double
epsilonr(double T)
{
	return 5321/T + 233.76 - 0.9297*T + 1.417*T*T/1000 - 0.8292*T*T*T/1000000;
}


PRIVATE INLINE double
bjerrum_length(double T)
{
	return 167100.052/(T*epsilonr(T));
}


PRIVATE INLINE double
ionic_strength(double rho)
{
  return rho;
}


PRIVATE INLINE double
kappa_inv(double rho, double T)
{
  return 8.1284/(sqrt(bjerrum_length(T)*ionic_strength(rho)));
}


PRIVATE INLINE double
kappa(double rho, double T)
{
  return sqrt(bjerrum_length(T)*ionic_strength(rho))/8.1284;
}


PRIVATE INLINE double
tau_ds(double T, double Helical_Rise)
{
  double bjerrum_length_inv, helical_rise_inv;

  bjerrum_length_inv = 1. / bjerrum_length(T);
  helical_rise_inv = 1. / Helical_Rise;
  return MIN2(helical_rise_inv, bjerrum_length_inv);
}


PRIVATE INLINE double
tau_ss(double T, double backbonelen)
{
  double bjerrum_length_inv;

  bjerrum_length_inv = 1/bjerrum_length(T);
  return MIN2(1/backbonelen, bjerrum_length_inv);
}


PRIVATE INLINE double
pairing_salt_const(double T, double Helical_Rise)
{
  return 2*(GASCONST / 1000.)*T*bjerrum_length(T)*Helical_Rise*tau_ds(T, Helical_Rise)*tau_ds(T, Helical_Rise);
}


PRIVATE double
approx_hyper(double y)
{
	double a, b, c;

	a = 1/(pow(y,6.)/pow(2*PI,6.)+1);
	b = pow(y, 4.)/(36*pow(PI, 4.)) - pow(y, 3.)/(24*PI*PI) + y*y/(2*PI*PI) - y/2;
	c = log(2*PI/y) - 1.96351;

	return a*b + (1-a)*c;
}


PRIVATE double
loop_salt_aux(double kmlss, int L, double T, double backbonelen)
{
	double a, b;

	a = (GASCONST / 1000.) * T * bjerrum_length(T) * L * backbonelen * tau_ss(T, backbonelen) * tau_ss(T, backbonelen);
	b = log(kmlss) - log(PI/2) + Eular_const + approx_hyper(kmlss) + 1/kmlss * (1 - exp(-kmlss) + kmlss*expn(1, kmlss));

	return a*b*100;
}

PUBLIC double
vrna_salt_loop(int L, double rho, double T, double backbonelen)
{
	if (L == 0)
		return 0;
	double correction, kmlss, kmlss_ref;
	
	kmlss_ref = kappa(VRNA_MODEL_DEFAULT_SALT, T) * L * backbonelen;
	kmlss = kappa(rho, T) * L * backbonelen;

	correction = loop_salt_aux(kmlss, L, T, backbonelen) - loop_salt_aux(kmlss_ref, L, T, backbonelen);

	return correction;
}

PUBLIC int
vrna_salt_loop_int(int L, double rho, double T, double backbonelen)
{
	double correction;

	correction = vrna_salt_loop(L, rho, T, backbonelen);
	return roundint(correction);
}



PUBLIC int
vrna_salt_stack(double rho, double T, double hrise)
{
	double correction, kn_ref, kappa_ref;
	
	kappa_ref = kappa(VRNA_MODEL_DEFAULT_SALT, T);
	kn_ref = kn(0, Rods_dist*kappa_ref);
	correction = 100*pairing_salt_const(T, hrise)*(kn(0, Rods_dist*kappa(rho, T)) - kn_ref);
	return roundint(correction);
}

/* Adapted from https://stackoverflow.com/a/19040841 */
/* LoopEnergy = m*size + b */
PUBLIC void
vrna_salt_ml(double saltLoop[], int lower, int upper, int* m, int* b)
{
	int sumx, sumxx;
	double y, sumy, sumyy, sumxy, denom, dm, db;

	sumx = sumxx = 0;
	sumy = sumyy = sumxy = 0.;

	for (int i=lower; i<=upper; i++)
	{ 
		sumx += i;
		sumxx += i * i;
		
		y = saltLoop[i];

		sumxy += i * y;
		sumy  += y;    
		sumyy += y*y;
	}

	denom = (double) ((upper-lower+1)*sumxx - sumx*sumx);
	dm = ((upper-lower+1) * sumxy  -  sumx * sumy) / denom;
	db = (sumy * sumxx  -  sumx * sumxy) / denom;

	*m = roundint(dm);
	*b = roundint(db);
	
}


/* Function to compute the salt correction for duplex initialization */
/* Fitted from 18 duplexes data (Chen & Znosko, 2013) */
PUBLIC int
vrna_salt_duplex_init(vrna_md_t *md_p)
{
  double    x;
  vrna_md_t md;

  if (md_p == NULL) {
    vrna_md_set_default(&md);
    md_p = &md;
  }

  if (md_p->saltDPXInit != 99999) {
    return md_p->saltDPXInit;
  } else {
    x = log(md_p->salt / VRNA_MODEL_DEFAULT_SALT);
    /* Converged duplex init correction */
    /* a = -1.25480589; */
    /* double b = -0.05306256; */
    /* int c = 160; */
    /* penalty = -exp(a*x*x+b*x+log(c)) + c; */
    /* a = -100.14040781; */
    return roundint(x * md_p->saltDPXInitFact);
  }
}
