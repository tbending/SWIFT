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

#ifndef SWIFT_LIGHTCONE_SHELL_H
#define SWIFT_LIGHTCONE_SHELL_H

/* Standard headers */
#include <stdio.h>

/* Config parameters. */
#include "../config.h"

/* Local headers */
#include "cosmology.h"
#include "lightcone_map.h"
#include "lightcone_map_types.h"
#include "particle_buffer.h"


enum lightcone_shell_state {
  shell_uninitialized,
  shell_current,
  shell_complete,
};


/**
 * @brief Information about a particle type contributing to the lightcone
 */
struct lightcone_particle_type {
  
  /*! Number of lightcone maps this particle type contributes to */
  int nr_maps;
  
  /*! Indices of the lightcone maps this particle type contributes to */
  int *map_index;

  /*! Amount of data to store per particle: theta, phi, radius and the value to add to each healpix map */
  size_t buffer_element_size;

};


/**
 * @brief Information about each lightcone shell
 */
struct lightcone_shell {

  /*! State of this shell */
  enum lightcone_shell_state state;

  /*! Inner radius of shell */
  double rmin;

  /*! Outer radius of shell */
  double rmax;

  /*! Minimum expansion factor for this shell */
  double amin;

  /*! Maximum expansion factor for this shell */
  double amax;

  /*! Number of maps associated with this shell */
  int nr_maps;

  /*! Array of lightcone maps for this shell */
  struct lightcone_map *map;

  /*! Buffers to store the map updates for each particle type */
  struct particle_buffer buffer[swift_type_count];

};


struct lightcone_shell *lightcone_shell_array_init(const struct cosmology *cosmo,
                                                   const char *radius_file, int nr_maps,
                                                   struct lightcone_map_type *map_type,
                                                   int nside, size_t total_nr_pix,
                                                   struct lightcone_particle_type *part_type,
                                                   size_t elements_per_block,
                                                   int *nr_shells_out);

void lightcone_shell_array_free(struct lightcone_shell *shell, int nr_shells);

void lightcone_shell_array_dump(const struct lightcone_shell *shell, int nr_shells, FILE *stream);

struct lightcone_shell *lightcone_shell_array_restore(FILE *stream, int nr_shells,
                                                      struct lightcone_particle_type *part_type,
                                                      size_t elements_per_block);

#endif /* SWIFT_LIGHTCONE_SHELL_H */
