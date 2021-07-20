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

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "black_holes.h"
#include "engine.h"
#include "gravity.h"
#include "hydro.h"
#include "lightcone_map.h"
#include "part.h"
#include "stars.h"

/* This object's header */
#include "lightcone_map_types.h"


/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_total_mass_type_contributes(int ptype) {

  switch(ptype) {
  case swift_type_gas:
  case swift_type_stars:
  case swift_type_black_hole:
  case swift_type_dark_matter:
  case swift_type_dark_matter_background:
  case swift_type_neutrino:
    return 1;
  default:
    return 0;
  }
}

/**
 * @brief Make a healpix map of projected mass in each pixel
 *
 * @param e the #engine structure
 * @param lightcone_props properties of the lightcone to update
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_total_mass_get_value(const struct engine *e,
                                          const struct lightcone_props *lightcone_props,
                                          const struct gpart *gp, const double a_cross,
                                          const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  /* const struct xpart *xparts = s->xparts; */ /* Currently not used */
  const struct spart *sparts = s->sparts;
  const struct bpart *bparts = s->bparts;

  switch (gp->type) {
  case swift_type_gas: {
    const struct part *p = &parts[-gp->id_or_neg_offset];
    return p->mass;
  } break;
  case swift_type_stars: {
    const struct spart *sp = &sparts[-gp->id_or_neg_offset];
    return sp->mass;
  } break;
  case swift_type_black_hole: {      
    const struct bpart *bp = &bparts[-gp->id_or_neg_offset];
    return bp->mass;
  } break;
  case swift_type_dark_matter:
  case swift_type_dark_matter_background:
  case swift_type_neutrino: {
    return gp->mass;
  } break;
  default:
    error("lightcone map function called on wrong particle type");
    return -1.0; /* Prevent 'missing return' error */
  }
}

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_gas_mass_type_contributes(int ptype) {

  switch(ptype) {
  case swift_type_gas:
    return 1;
  default:
    return 0;
  }
}

/**
 * @brief Make a healpix map of projected gas mass in each pixel
 *
 * @param e the #engine structure
 * @param lightcone_props properties of the lightcone to update
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_gas_mass_get_value(const struct engine *e,
                                        const struct lightcone_props *lightcone_props,
                                        const struct gpart *gp, const double a_cross,
                                        const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;

  switch (gp->type) {
  case swift_type_gas: {
    const struct part *p = &parts[-gp->id_or_neg_offset];
    return p->mass;
  } break;
  default:
    error("lightcone map function called on wrong particle type");
    return -1.0;  /* Prevent 'missing return' error */
  }
}

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_dark_matter_mass_type_contributes(int ptype) {

  switch(ptype) {
  case swift_type_dark_matter:
  case swift_type_dark_matter_background:
    return 1;
  default:
    return 0;
  }
}

/**
 * @brief Make a healpix map of projected dark matter mass in each pixel
 *
 * @param e the #engine structure
 * @param lightcone_props properties of the lightcone to update
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_dark_matter_mass_get_value(const struct engine *e,
                                                const struct lightcone_props *lightcone_props,
                                                const struct gpart *gp, const double a_cross,
                                                const double x_cross[3]) {
  switch (gp->type) {
  case swift_type_dark_matter:
  case swift_type_dark_matter_background: {
    return gp->mass;
  } break;
  default:
    error("lightcone map function called on wrong particle type");
    return -1.0;  /* Prevent 'missing return' error */
  }
}

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_stellar_mass_type_contributes(int ptype) {

  switch(ptype) {
  case swift_type_stars:
    return 1;
  default:
    return 0;
  }
}

/**
 * @brief Make a healpix map of stellar mass in each pixel
 *
 * @param e the #engine structure
 * @param lightcone_props properties of the lightcone to update
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_stellar_mass_get_value(const struct engine *e,
                                            const struct lightcone_props *lightcone_props,
                                            const struct gpart *gp, const double a_cross,
                                            const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct spart *sparts = s->sparts;

  switch (gp->type) {
  case swift_type_stars: {
    const struct spart *sp = &sparts[-gp->id_or_neg_offset];
    return sp->mass;
  } break;
  default:
    error("lightcone map function called on wrong particle type");
    return -1.0; /* Prevent 'missing return' error */
  }
}

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_neutrino_mass_type_contributes(int ptype) {

  switch(ptype) {
  case swift_type_neutrino:
    return 1;
  default:
    return 0;
  }
}

/**
 * @brief Make a healpix map of projected neutrino mass in each pixel
 *
 * @param e the #engine structure
 * @param lightcone_props properties of the lightcone to update
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_neutrino_mass_get_value(const struct engine *e,
                                             const struct lightcone_props *lightcone_props,
                                             const struct gpart *gp, const double a_cross,
                                             const double x_cross[3]) {
  switch (gp->type) {
  case swift_type_neutrino: {
    return gp->mass;
  } break;
  default:
    error("lightcone map function called on wrong particle type");
    return -1.0;  /* Prevent 'missing return' error */
  }
}

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_compton_y_type_contributes(int ptype) {

  switch(ptype) {
  case swift_type_gas:
    return 1;
  default:
    return 0;
  }
}

/**
 * @brief Make a healpix map of the compton y parameter
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_compton_y_get_value(const struct engine *e,
                                         const struct lightcone_props *lightcone_props,
                                         const struct gpart *gp, const double a_cross,
                                         const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart* xparts = e->s->xparts;

  /* Handle on the physics modules */
  const struct cosmology* cosmo = e->cosmology;
  const struct hydro_props* hydro_props = e->hydro_properties;
  const struct unit_system* us = e->internal_units;
  const struct phys_const* phys_const = e->physical_constants;
  const struct cooling_function_data* cool_func = e->cooling_func;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart* xp = &xparts[-gp->id_or_neg_offset];
      double y_compton = cooling_get_ycompton(phys_const, hydro_props, us, cosmo,
                                       cool_func, p, xp);
      double x_squared = x_cross[0] * x_cross[0] * a_cross * a_cross;
      double y_squared = x_cross[1] * x_cross[1] * a_cross * a_cross;
      double z_squared = x_cross[2] * x_cross[2] * a_cross * a_cross;
      double angular_diameter_distance_2 = x_squared + y_squared + z_squared;

      double pixel_size_2 = lightcone_props->pixel_area_steradians;

      double y_for_map = y_compton / (pixel_size_2 * angular_diameter_distance_2);

      return y_for_map;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0;  /* Prevent 'missing return' error */
      break;
  }
}

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_doppler_b_type_contributes(int ptype) {

  switch(ptype) {
  case swift_type_gas:
    return 1;
  default:
    return 0;
  }
}

/**
 * @brief Make a healpix map of the doppler b parameter
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_doppler_b_get_value(const struct engine *e,
                                         const struct lightcone_props *lightcone_props,
                                         const struct gpart *gp, const double a_cross,
                                         const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart* xparts = e->s->xparts;

  /* Handle on the physics modules */
  const struct cosmology* cosmo = e->cosmology;
  const struct hydro_props* hydro_props = e->hydro_properties;
  const struct unit_system* us = e->internal_units;
  const struct phys_const* phys_const = e->physical_constants;
  const struct cooling_function_data* cool_func = e->cooling_func;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart* xp = &xparts[-gp->id_or_neg_offset];
      double n_e = cooling_get_electron_pressure(phys_const, hydro_props, us, cosmo,
                                          cool_func, p, xp);

      double rho = hydro_get_physical_density(p, cosmo);

      double m = hydro_get_mass(p);

      const double c = phys_const->const_speed_light_c;

      const double sigma_thompson = phys_const->const_thomson_cross_section;

      double doppler_b_factor = n_e * m * sigma_thompson / (rho * c);

      double x_squared = x_cross[0] * x_cross[0] * a_cross * a_cross;
      double y_squared = x_cross[1] * x_cross[1] * a_cross * a_cross;
      double z_squared = x_cross[2] * x_cross[2] * a_cross * a_cross;
      double angular_diameter_distance_2 = x_squared + y_squared + z_squared;

      double angular_diameter_distance = sqrt(angular_diameter_distance_2);

      double radial_velocity =
          (p->v[0] * x_cross[0] * a_cross + p->v[1] * x_cross[1] * a_cross +
           p->v[2] * x_cross[2] * a_cross) /
          angular_diameter_distance;

      double pixel_size_2 = lightcone_props->pixel_area_steradians;

      double b_for_map = doppler_b_factor * radial_velocity /
                  (pixel_size_2 * angular_diameter_distance_2);

      return b_for_map;
    } break;
    default:
      error("lightcone map function called on wrong particle type");
      return -1.0;  /* Prevent 'missing return' error */
      break;
  }
}

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_dispersion_meassure_type_contributes(int ptype) {

  switch(ptype) {
  case swift_type_gas:
    return 1;
  default:
    return 0;
  }
}

/**
 * @brief Make a healpix map of the dispersion meassure
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the
 * lightcone
 */
double lightcone_map_dispersion_meassure_get_value(const struct engine *e,
                                                   const struct lightcone_props *lightcone_props,
                                                   const struct gpart *gp, const double a_cross,
                                                   const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart* xparts = e->s->xparts;

  /* Handle on the physics modules */
  const struct cosmology* cosmo = e->cosmology;
  const struct hydro_props* hydro_props = e->hydro_properties;
  const struct unit_system* us = e->internal_units;
  const struct phys_const* phys_const = e->physical_constants;
  const struct cooling_function_data* cool_func = e->cooling_func;

  switch (gp->type) {
    case swift_type_gas: {
      const struct part *p = &parts[-gp->id_or_neg_offset];
      const struct xpart* xp = &xparts[-gp->id_or_neg_offset];
      double n_e = cooling_get_electron_pressure(phys_const, hydro_props, us, cosmo,
                                          cool_func, p, xp);

      double rho = hydro_get_physical_density(p, cosmo);

      double m = hydro_get_mass(p);

      double number_of_electrons = n_e * m / (rho);

      double x_squared = x_cross[0] * x_cross[0] * a_cross * a_cross;
      double y_squared = x_cross[1] * x_cross[1] * a_cross * a_cross;
      double z_squared = x_cross[2] * x_cross[2] * a_cross * a_cross;
      double angular_diameter_distance_2 = x_squared + y_squared + z_squared;

      double pixel_size_2 = lightcone_props->pixel_area_steradians;

      double dm_for_map = number_of_electrons /
                  (pixel_size_2 * angular_diameter_distance_2);

      return dm_for_map;
    } break;
    default:
      /* Not gas, nothing to do */
      error("lightcone map function called on wrong particle type");
      return -1.0;  /* Prevent 'missing return' error */
      break;
  }
}
