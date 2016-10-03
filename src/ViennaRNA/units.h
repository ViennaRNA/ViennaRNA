#ifndef VIENNA_RNA_PACKAGE_UNITS_H
#define VIENNA_RNA_PACKAGE_UNITS_H

/**
 *  @file     units.h
 *  @ingroup  units
 *  @brief    Physical Units and Functions to convert them into each other
 */

/**
 *  @{
 *  @ingroup  units
 *
 *
 */

/**
 *  @brief Energy / Work Units
 *
 *  @see vrna_convert_energy()
 */
typedef enum {
  VRNA_UNIT_J,            /**< @brief Joule ( @f$ 1~J = 1~kg \cdot m^2 s^{-2} @f$ ) */
  VRNA_UNIT_KJ,           /**< @brief Kilojoule ( @f$ 1~kJ = 1,000~J @f$ ) */
  VRNA_UNIT_CAL_IT,       /**< @brief Calorie (International (Steam) Table, @f$ 1~cal_{IT} = 4.1868~J @f$ ) */
  VRNA_UNIT_DACAL_IT,     /**< @brief Decacolorie (International (Steam) Table, @f$ 1~dacal_{IT} = 10~cal_{IT} = 41.868~J @f$ ) */
  VRNA_UNIT_KCAL_IT,      /**< @brief Kilocalorie (International (Steam) Table, @f$ 1~kcal_{IT} = 4.1868~kJ @f$ ) */
  VRNA_UNIT_CAL,          /**< @brief Calorie (Thermochemical, @f$ 1~cal_{th} = 4.184~J @f$ ) */
  VRNA_UNIT_DACAL,        /**< @brief Decacalorie (Thermochemical, @f$ 1~dacal_{th} = 10~cal_{th} = 41.84~J @f$ ) */
  VRNA_UNIT_KCAL,         /**< @brief Kilocalorie (Thermochemical, @f$ 1~kcal_{th} = 4.184~kJ @f$ ) */
  VRNA_UNIT_G_TNT,        /**< @brief g TNT ( @f$ 1~g~\mathrm{TNT} = 1,000~cal_{th} = 4,184~J @f$ ) */
  VRNA_UNIT_KG_TNT,       /**< @brief kg TNT ( @f$ 1~kg~\mathrm{TNT} = 1,000~kcal_{th} = 4,184~kJ @f$ ) */
  VRNA_UNIT_T_TNT,        /**< @brief ton TNT ( @f$ 1~t~\mathrm{TNT} = 1,000,000~kcal_{th} = 4,184~MJ @f$ ) */
  VRNA_UNIT_EV,           /**< @brief Electronvolt ( @f$ 1~eV = 1.602176565 \times 10^{-19}~J @f$ ) */
  VRNA_UNIT_WH,           /**< @brief Watt hour ( @f$ 1~W \cdot h = 1~W \cdot 3,600 s = 3,600~J = 3.6~kJ @f$ ) */
  VRNA_UNIT_KWH,          /**< @brief Kilowatt hour ( @f$ 1~kW \cdot h = 1~kW \cdot 3,600~s = 3,600~kJ = 3.6~MJ @f$ ) */
} vrna_unit_energy_e;


/**
 *  @brief Temperature Units
 *
 *  @see vrna_convert_temperature()
 */
typedef enum {
  VRNA_UNIT_K,            /**< @brief Kelvin (K) */
  VRNA_UNIT_DEG_C,        /**< @brief Degree Celcius (&deg;C) (@f$ [^{\circ}C] = [K] - 273.15 @f$ ) */
  VRNA_UNIT_DEG_F,        /**< @brief Degree Fahrenheit (&deg;F) (@f$ [^{\circ}F] = [K] \times \frac{9}{5} - 459.67 @f$ ) */
  VRNA_UNIT_DEG_R,        /**< @brief Degree Rankine (&deg;R) (@f$ [^{\circ}R] = [K] \times \frac{9}{5} @f$ ) */
  VRNA_UNIT_DEG_N,        /**< @brief Degree Newton (&deg;N) (@f$ [^{\circ}N] = ([K] - 273.15) \times \frac{33}{100} @f$ ) */
  VRNA_UNIT_DEG_DE,       /**< @brief Degree Delisle (&deg;De) (@f$ [^{\circ}De] = (373.15 - [K]) \time \frac{3}{2} @f$ ) */
  VRNA_UNIT_DEG_RE,       /**< @brief Degree R&eacute;aumur (&deg;R&eacute;) (@f$ [^{\circ}R{\acute e}] = ([K] - 273.15) \times \frac{4}{5} @f$ ) */
  VRNA_UNIT_DEG_RO,       /**< @brief Degree R&oslash;mer (&deg;R&oslash;) (@f$ [^{\circ}\mathrm{R{\o}}] = ([K] - 273.15) \times \frac{21}{40} + 7.5 @f$ ) */
} vrna_unit_temperature_e;


/**
 *  @brief  Convert between energy / work units
 *
 *  @see #vrna_unit_energy_e
 *  @param  energy  Input energy value
 *  @param  from    Input unit
 *  @param  to      Output unit
 *  @return         Energy value in Output unit
 */
double vrna_convert_energy( double energy,
                            vrna_unit_energy_e from,
                            vrna_unit_energy_e to);


/**
 *  @brief  Convert between temperature units
 *
 *  @see #vrna_unit_temperature_e
 *  @param  temp    Input temperature value
 *  @param  from    Input unit
 *  @param  to      Output unit
 *  @return         Temperature value in Output unit
 */
double vrna_convert_temperature(double temp,
                                vrna_unit_temperature_e from,
                                vrna_unit_temperature_e to);

/**
 *  @}
 */

#endif
