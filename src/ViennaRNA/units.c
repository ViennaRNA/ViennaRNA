/*
            units.c

   (c) 2016 Ronny Lorenz
            ViennaRNA Package
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include "ViennaRNA/utils.h"
#include "ViennaRNA/units.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/*
#################################
# GLOBAL VARIABLES              #
#################################
*/


/*
#################################
# PRIVATE VARIABLES             #
#################################
*/


/*
#################################
# PRIVATE FUNCTION DECLARATIONS #
#################################
*/
PRIVATE INLINE double Joule2CalIT(double e);
PRIVATE INLINE double CalIT2Joule(double e);
PRIVATE INLINE double Joule2Cal(double e);
PRIVATE INLINE double Cal2Joule(double e);
PRIVATE INLINE double Joule2TNT(double e);
PRIVATE INLINE double TNT2Joule(double e);
PRIVATE INLINE double ElectronVolt2Joule(double e);
PRIVATE INLINE double Joule2ElectronVolt(double e);
PRIVATE INLINE double Watthour2Joule(double e);
PRIVATE INLINE double Joule2Watthour(double e);

/*
#################################
# BEGIN OF FUNCTION DEFINITIONS #
#################################
*/

PUBLIC double
vrna_convert_energy(double input,
                    vrna_unit_energy_e from,
                    vrna_unit_energy_e to){

  double output;

  /* we first convert the input to 1 kJ and then convert to the target unit */

  switch(from){
    /* International (Steam) Table where 1 cal_IT = 4.1868 J */
    case VRNA_UNIT_CAL_IT:      input /= 10.;
                                /* fall through */
    case VRNA_UNIT_DACAL_IT:    input /= 100.;
                                /* fall through */
    case VRNA_UNIT_KCAL_IT:     input = CalIT2Joule(input);
                                break;

    /* Thermochemicalrie 1 cal_TH = 4.184 J */
    case VRNA_UNIT_CAL:         input /= 10.;
                                /* fall through */
    case VRNA_UNIT_DACAL:       input /= 100.;
                                /* fall through */
    case VRNA_UNIT_KCAL:        input = Cal2Joule(input);
                                break;

    case VRNA_UNIT_J:           input /= 1000.;
    case VRNA_UNIT_KJ:          break;

    /* kg TNT equivalent 1 kg TNT = 1,000 kcal_{th} = 4,184 kJ */
    case VRNA_UNIT_G_TNT:       input /= 1000.;
                                /* fall through */
    case VRNA_UNIT_KG_TNT:      input /= 1000.;
                                /* fall through */
    case VRNA_UNIT_T_TNT:       input = TNT2Joule(input);
                                break;

    case VRNA_UNIT_EV:          input = ElectronVolt2Joule(input);
                                break;

    case VRNA_UNIT_WH:          input /= 1000.;
                                /* fall through */
    case VRNA_UNIT_KWH:         input = Watthour2Joule(input);
                                break;

    /* we assume kcal_TH by default */
    default:                    input = Cal2Joule(input);
                                break;
  }

  /* Lets convert to target unit now */
  switch(to){
    case VRNA_UNIT_J:           input *= 1000.;
                                /* fall through */
    case VRNA_UNIT_KJ:          output = input;
                                break;

    case VRNA_UNIT_CAL:         input *= 10.;
                                /* fall through */
    case VRNA_UNIT_DACAL:       input *= 100.;
                                /* fall through */
    case VRNA_UNIT_KCAL:        output = Joule2Cal(input);
                                break;

    case VRNA_UNIT_CAL_IT:      input *= 10.;
                                /* fall through */
    case VRNA_UNIT_DACAL_IT:    input *= 100.;
                                /* fall through */
    case VRNA_UNIT_KCAL_IT:     output = Joule2CalIT(input);
                                break;

    case VRNA_UNIT_G_TNT:       input *= 1000.;
                                /* fall through */
    case VRNA_UNIT_KG_TNT:      input *= 1000.;
                                /* fall through */
    case VRNA_UNIT_T_TNT:       output = Joule2TNT(input);
                                break;

    case VRNA_UNIT_EV:          output = Joule2ElectronVolt(input);
                                break;

    case VRNA_UNIT_WH:          input *= 1000.;
                                /* fall through */
    case VRNA_UNIT_KWH:         output = Joule2Watthour(input);
                                break;

    /* we assume kcal_TH by default */
    default:                    output = Joule2Cal(input);
                                break;
  }

  return output;
}


PUBLIC double
vrna_convert_temperature( double temp,
                          vrna_unit_temperature_e from,
                          vrna_unit_temperature_e to){

  switch(from){
    case VRNA_UNIT_DEG_C:   temp += K0;
                            break;
    case VRNA_UNIT_DEG_F:   temp += 459.67;
                            temp *= 5. / 9.;
                            break;
    case VRNA_UNIT_DEG_R:   temp /= 9. / 5.;
                            break;
    case VRNA_UNIT_DEG_N:   temp *= 100. / 33.;
                            temp += K0;
                            break;
    case VRNA_UNIT_DEG_DE:  temp *= 2. / 3.;
                            temp = 100. + K0 - temp;
                            break;
    case VRNA_UNIT_DEG_RE:  temp *= 5. / 4.;
                            temp += K0;
                            break;
    case VRNA_UNIT_DEG_RO:  temp -= 7.5;
                            temp *= 40. / 21.;
                            temp += K0;
                            break;
    case VRNA_UNIT_K:       /* fall through */
    default:                break;
  }

  switch(to){
    case VRNA_UNIT_DEG_C:   temp -= K0;
                            break;
    case VRNA_UNIT_DEG_F:   temp *= 9. / 5.;
                            temp -= 459.67;
                            break;
    case VRNA_UNIT_DEG_R:   temp *= 9. / 5.;
                            break;
    case VRNA_UNIT_DEG_N:   temp -= K0;
                            temp *= 33. / 100.;
                            break;
    case VRNA_UNIT_DEG_DE:  temp = 100. + K0 - temp;
                            temp *= 3. / 2.;
                            break;
    case VRNA_UNIT_DEG_RE:  temp -= K0;
                            temp *= 4. / 5.;
                            break;
    case VRNA_UNIT_DEG_RO:  temp -= K0;
                            temp *= 21. / 40.;
                            temp += 7.5;
                            break;
    case VRNA_UNIT_K:       /* fall through */
    default:                break;
  }

  return temp;
}

/*
#################################
# STATIC helper functions below #
#################################
*/

PRIVATE INLINE double
Joule2CalIT(double e){

  return (e * 4.1868);
}


PRIVATE INLINE double
CalIT2Joule(double e){

  return (e / 4.1868);
}


PRIVATE INLINE double
Joule2Cal(double e){

  return (e * 4.184);
}


PRIVATE INLINE double
Cal2Joule(double e){

  return (e / 4.184);
}


PRIVATE INLINE double
Joule2TNT(double e){

  /* 1 kg TNT = 1000 kcal_TH = 4184 KJ */
  return (e / 4184000.);
}


PRIVATE INLINE double
TNT2Joule(double e){

  /* 1 kg TNT = 1000 kcal_TH = 4184 KJ */
  return (e * 4184000.);
}


PRIVATE INLINE double
ElectronVolt2Joule(double e){

  return (1.602176565e-19 * e);
}


PRIVATE INLINE double
Joule2ElectronVolt(double e){

  return (6.24150974e18 * e);
}


PRIVATE INLINE double
Watthour2Joule(double e){

  return (3600 * e);
}


PRIVATE INLINE double
Joule2Watthour(double e){

  return (e / 3600);
}
