/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 John Helly (j.c.helly@durham.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef SWIFT_LIGHTCONE_MAP_TYPES_H
#define SWIFT_LIGHTCONE_MAP_TYPES_H

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "parser.h"
#include "part_type.h"
#include "units.h"

/* Avoid cyclic inclusions */
struct engine;
struct lightcone_map;
struct gpart;

/* Type to store pointer to function for updating a healpix map */
typedef double (*map_update_function_t)(const struct engine *e,
                                        const struct gpart *gp, const double a_cross,
                                        const double x_cross[3]);

/* Type to store pointer to function to check which types contribute to a map */
typedef int (*map_contrib_function_t)(int ptype);

/**
 * @brief Struct to store information on one type of lightcone map
 */
struct lightcone_map_type {
  char name[PARSER_MAX_LINE_SIZE];
  map_update_function_t update_map;
  map_contrib_function_t ptype_contributes;
  enum unit_conversion_factor units;
};


/* 
   Healpix map of total mass
*/
int lightcone_map_total_mass_type_contributes(int ptype);

double lightcone_map_total_mass_get_value(const struct engine *e,
                                          const struct gpart *gp, const double a_cross,
                                          const double x_cross[3]);
/* 
   Healpix map of gas mass
*/
int lightcone_map_gas_mass_type_contributes(int ptype);

double lightcone_map_gas_mass_get_value(const struct engine *e,
                                        const struct gpart *gp, const double a_cross,
                                        const double x_cross[3]);
/* 
   Healpix map of dark matter mass
*/
int lightcone_map_dark_matter_mass_type_contributes(int ptype);

double lightcone_map_dark_matter_mass_get_value(const struct engine *e,
                                                const struct gpart *gp, const double a_cross,
                                                const double x_cross[3]);
/* 
   Healpix map of stellar mass
*/
int lightcone_map_stellar_mass_type_contributes(int ptype);

double lightcone_map_stellar_mass_get_value(const struct engine *e,
                                                const struct gpart *gp, const double a_cross,
                                                const double x_cross[3]);
/* 
   Healpix map of neutrino mass
*/
int lightcone_map_neutrino_mass_type_contributes(int ptype);

double lightcone_map_neutrino_mass_get_value(const struct engine *e,
                                             const struct gpart *gp, const double a_cross,
                                             const double x_cross[3]);
/* 
   Healpix map of compton y
*/
int lightcone_map_compton_y_type_contributes(int ptype);

double lightcone_map_compton_y_get_value(const struct engine *e,
                                       const struct gpart *gp, const double a_cross,
                                       const double x_cross[3]);
/* 
   Healpix map of doppler b
*/
int lightcone_map_doppler_b_type_contributes(int ptype);

double lightcone_map_doppler_b_get_value(const struct engine *e,
                                       const struct gpart *gp, const double a_cross,
                                       const double x_cross[3]);
/* 
   Healpix map of dispersion meassure
*/
int lightcone_map_dispersion_meassure_type_contributes(int ptype);

double lightcone_map_dispersion_meassure_get_value(const struct engine *e,
                                                 const struct gpart *gp, const double a_cross,
                                                 const double x_cross[3]);


/* This associates map names to the appropriate update function and unit info */
static const struct lightcone_map_type lightcone_map_types[] = {
  {"TotalMass",      lightcone_map_total_mass_get_value,       lightcone_map_total_mass_type_contributes,       UNIT_CONV_MASS},
  {"GasMass",        lightcone_map_gas_mass_get_value,         lightcone_map_gas_mass_type_contributes,         UNIT_CONV_MASS},
  {"DarkMatterMass", lightcone_map_dark_matter_mass_get_value, lightcone_map_dark_matter_mass_type_contributes, UNIT_CONV_MASS},
  {"StellarMass",    lightcone_map_stellar_mass_get_value,     lightcone_map_stellar_mass_type_contributes,     UNIT_CONV_MASS},
  {"NeutrinoMass",   lightcone_map_neutrino_mass_get_value,    lightcone_map_neutrino_mass_type_contributes,    UNIT_CONV_MASS},
  {"ComptonY",       lightcone_map_compton_y_get_value,        lightcone_map_compton_y_type_contributes,    UNIT_CONV_NO_UNITS},
  {"DopplerB",       lightcone_map_doppler_b_get_value,        lightcone_map_doppler_b_type_contributes,    UNIT_CONV_NO_UNITS},
  {"DM",             lightcone_map_dispersion_meassure_get_value, lightcone_map_dispersion_meassure_type_contributes, UNIT_CONV_INV_AREA},
  {"",               NULL,                                     NULL,                                            UNIT_CONV_NO_UNITS},
  /* NULL functions indicate end of array */
};

#endif
