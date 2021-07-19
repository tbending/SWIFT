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


/* Required for the xrays */
#include "extra_io/EAGLE/extra.h"
#include "io_properties.h"


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
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_total_mass_get_value(const struct engine *e,
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
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_gas_mass_get_value(const struct engine *e,
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
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_dark_matter_mass_get_value(const struct engine *e,
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
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_stellar_mass_get_value(const struct engine *e,
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
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_neutrino_mass_get_value(const struct engine *e,
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


#ifndef SWIFT_EXTRA_IO_EAGLE_H
#define SWIFT_EXTRA_IO_EAGLE_H

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_xray_erosita_low_intrinsic_photons_type_contributes(int ptype) {

  switch(ptype) {
  case swift_type_gas:
    return 1;
  default:
    return 0;
  }
}


/**
 * @brief Make a healpix map of projected erosita-low intrinsic photon flux in each pixel
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_xray_erosita_low_intrinsic_photons_get_value(const struct engine *e,
                                        const struct gpart *gp, const double a_cross,
                                        const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;

  switch (gp->type) {
  case swift_type_gas: {
    const struct part *p = &parts[-gp->id_or_neg_offset];
    const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

    const double z_cross = (1 / a_cross) - 1;
    const double cdist_cross = sqrt( pow(x_cross[0], 2) + pow(x_cross[1], 2) + pow(x_cross[2], 2) );

    const double luminosity = extra_io_get_xray_fluxes(
    p, xp, e, xray_band_types_erosita_low_intrinsic_photons);

    const double flux = luminosity / (4 * M_PI * pow(cdist_cross, 2) * (1 + z_cross)); // photon luminosity distance

    return flux;
  } break;
  default:
    error("lightcone map function called on wrong particle type");
    return -1.0;  /* Prevent 'missing return' error */
  }
}


/* erosita-low energy flux */

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_xray_erosita_low_intrinsic_type_contributes(int ptype) {

  switch(ptype) {
  case swift_type_gas:
    return 1;
  default:
    return 0;
  }
}


/**
 * @brief Make a healpix map of projected erosita-low intrinsic energy flux in each pixel
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_xray_erosita_low_intrinsic_energy_get_value(const struct engine *e,
                                        const struct gpart *gp, const double a_cross,
                                        const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;

  switch (gp->type) {
  case swift_type_gas: {
    const struct part *p = &parts[-gp->id_or_neg_offset];
    const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

    const double z_cross = (1 / a_cross) - 1;
    const double cdist_cross = sqrt( pow(x_cross[0], 2) + pow(x_cross[1], 2) + pow(x_cross[2], 2) );

    const double luminosity = extra_io_get_xray_fluxes(
    p, xp, e, xray_band_types_erosita_low_intrinsic_energies);

    const double flux = luminosity / (4 * M_PI * pow(cdist_cross, 2) * pow((1 + z_cross), 2)); // energy luminosity distance

    return flux;
  } break;
  default:
    error("lightcone map function called on wrong particle type");
    return -1.0;  /* Prevent 'missing return' error */
  }
}

/* erosita_high photon flux */

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_xray_erosita_high_intrinsic_photons_type_contributes(int ptype) {

  switch(ptype) {
  case swift_type_gas:
    return 1;
  default:
    return 0;
  }
}


/**
 * @brief Make a healpix map of projected erosita-high intrinsic photon flux in each pixel
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_xray_erosita_high_intrinsic_photons_get_value(const struct engine *e,
                                        const struct gpart *gp, const double a_cross,
                                        const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;

  switch (gp->type) {
  case swift_type_gas: {
    const struct part *p = &parts[-gp->id_or_neg_offset];
    const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

    const double z_cross = (1 / a_cross) - 1;
    const double cdist_cross = sqrt( pow(x_cross[0], 2) + pow(x_cross[1], 2) + pow(x_cross[2], 2) );

    const double luminosity = extra_io_get_xray_fluxes(
    p, xp, e, xray_band_types_erosita_high_intrinsic_photons);

    const double flux = luminosity / (4 * M_PI * pow(cdist_cross, 2) * (1 + z_cross)); // photon luminosity distance

    return flux;
  } break;
  default:
    error("lightcone map function called on wrong particle type");
    return -1.0;  /* Prevent 'missing return' error */
  }
}


/* erosita-high energy flux */

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_xray_erosita_high_intrinsic_type_contributes(int ptype) {

  switch(ptype) {
  case swift_type_gas:
    return 1;
  default:
    return 0;
  }
}


/**
 * @brief Make a healpix map of projected erosita-high intrinsic energy flux in each pixel
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_xray_erosita_high_intrinsic_energy_get_value(const struct engine *e,
                                        const struct gpart *gp, const double a_cross,
                                        const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;

  switch (gp->type) {
  case swift_type_gas: {
    const struct part *p = &parts[-gp->id_or_neg_offset];
    const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

    const double z_cross = (1 / a_cross) - 1;
    const double cdist_cross = sqrt( pow(x_cross[0], 2) + pow(x_cross[1], 2) + pow(x_cross[2], 2) );

    const double luminosity = extra_io_get_xray_fluxes(
    p, xp, e, xray_band_types_erosita_high_intrinsic_energies);

    const double flux = luminosity / (4 * M_PI * pow(cdist_cross, 2) * pow((1 + z_cross), 2)); // energy luminosity distance

    return flux;
  } break;
  default:
    error("lightcone map function called on wrong particle type");
    return -1.0;  /* Prevent 'missing return' error */
  }
}

/* ROSAT photon flux */

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_xray_rosat_intrinsic_photons_type_contributes(int ptype) {

  switch(ptype) {
  case swift_type_gas:
    return 1;
  default:
    return 0;
  }
}


/**
 * @brief Make a healpix map of projected ROSAT intrinsic photon flux in each pixel
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_xray_rosat_intrinsic_photons_get_value(const struct engine *e,
                                        const struct gpart *gp, const double a_cross,
                                        const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;

  switch (gp->type) {
  case swift_type_gas: {
    const struct part *p = &parts[-gp->id_or_neg_offset];
    const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

    const double z_cross = (1 / a_cross) - 1;
    const double cdist_cross = sqrt( pow(x_cross[0], 2) + pow(x_cross[1], 2) + pow(x_cross[2], 2) );

    const double luminosity = extra_io_get_xray_fluxes(
    p, xp, e, xray_band_types_ROSAT_intrinsic_photons);

    const double flux = luminosity / (4 * M_PI * pow(cdist_cross, 2) * (1 + z_cross)); // photon luminosity distance

    return flux;
  } break;
  default:
    error("lightcone map function called on wrong particle type");
    return -1.0;  /* Prevent 'missing return' error */
  }
}


/* ROSAT energy flux */

/**
 * @brief Determine if a particle type contributes to this map type
 *
 * @param part_type the particle type
 */
int lightcone_map_xray_rosat_intrinsic_type_contributes(int ptype) {

  switch(ptype) {
  case swift_type_gas:
    return 1;
  default:
    return 0;
  }
}


/**
 * @brief Make a healpix map of projected ROSAT intrinsic energy flux in each pixel
 *
 * @param e the #engine structure
 * @param gp the #gpart to add to the map
 * @param a_cross expansion factor at which the particle crosses the lightcone
 * @param x_cross comoving coordinates at which the particle crosses the lightcone
 */
double lightcone_map_xray_rosat_intrinsic_energy_get_value(const struct engine *e,
                                        const struct gpart *gp, const double a_cross,
                                        const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;

  switch (gp->type) {
  case swift_type_gas: {
    const struct part *p = &parts[-gp->id_or_neg_offset];
    const struct xpart *xp = &xparts[-gp->id_or_neg_offset];

    const double z_cross = (1 / a_cross) - 1;
    const double cdist_cross = sqrt( pow(x_cross[0], 2) + pow(x_cross[1], 2) + pow(x_cross[2], 2) );

    const double luminosity = extra_io_get_xray_fluxes(
    p, xp, e, xray_band_types_ROSAT_intrinsic_energies);

    const double flux = luminosity / (4 * M_PI * pow(cdist_cross, 2) * pow((1 + z_cross), 2)); // energy luminosity distance

    return flux;
  } break;
  default:
    error("lightcone map function called on wrong particle type");
    return -1.0;  /* Prevent 'missing return' error */
  }
}

#endif /* SWIFT_EXTRA_IO_EAGLE_H */

