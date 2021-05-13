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

/* Healpix C API */
#ifdef HAVE_CHEALPIX
#include <chealpix.h>
#endif


void lightcone_map_total_mass(struct lightcone_map *map, const struct engine *e,
                              const struct gpart *gp, const double a_cross,
                              const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  /* const struct xpart *xparts = s->xparts; */ /* Currently not used */
  const struct spart *sparts = s->sparts;
  const struct bpart *bparts = s->bparts;

  /* Angular smoothing radius for this particle */
  const double radius = 0.0;

  switch (gp->type) {
  case swift_type_gas: {
    const struct part *p = &parts[-gp->id_or_neg_offset];
    lightcone_map_buffer_update(map, x_cross, radius, p->mass);
  } break;
  case swift_type_stars: {
    const struct spart *sp = &sparts[-gp->id_or_neg_offset];
    lightcone_map_buffer_update(map, x_cross, radius, sp->mass);
  } break;
  case swift_type_black_hole: {      
    const struct bpart *bp = &bparts[-gp->id_or_neg_offset];
    lightcone_map_buffer_update(map, x_cross, radius, bp->mass);
  } break;
  case swift_type_dark_matter:
  case swift_type_dark_matter_background:
  case swift_type_neutrino: {
    lightcone_map_buffer_update(map, x_cross, radius, gp->mass);
  } break;
  default:
    /* Unknown type, nothing to do */
    break;
  }
}


void lightcone_map_gas_mass(struct lightcone_map *map, const struct engine *e,
                            const struct gpart *gp, const double a_cross,
                            const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;

  /* Angular smoothing radius for this particle */
  const double radius = 0.0;

  switch (gp->type) {
  case swift_type_gas: {
    const struct part *p = &parts[-gp->id_or_neg_offset];
    lightcone_map_buffer_update(map, x_cross, radius, p->mass);
  } break;
  default:
    /* Not gas, nothing to do */
    break;
  }
}


void lightcone_map_neutrino_mass(struct lightcone_map *map, const struct engine *e,
                                 const struct gpart *gp, const double a_cross,
                                 const double x_cross[3]) {

  /* Angular smoothing radius for this particle */
  const double radius = 0.0;
  
  switch (gp->type) {
  case swift_type_neutrino: {
    lightcone_map_buffer_update(map, x_cross, radius, gp->mass);
  } break;
  default:
    /* Not a neutrino, nothing to do */
    break;
  }
}

