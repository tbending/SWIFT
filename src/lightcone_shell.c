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

/* Some standard headers. */
#include <stdio.h>
#include <stdlib.h>

/* Local headers */
#include "cosmology.h"
#include "engine.h"

/* This object's header. */
#include "lightcone_shell.h"


/**
 * @brief Read in shell radii for lightcone healpix maps
 *
 * Allocates the output array, shell_out.
 *
 * @param cosmo the #cosmology structure
 * @param radius_file name of the file with shell radii
 * @param nr_shells returns the number of shells
 * @param shell_out returns the array of shells
 */
static void read_shell_radii(const struct cosmology *cosmo, const char *radius_file,
                             int *nr_shells, struct lightcone_shell **shell_out) {
  

  /* Allow shell radii to be specified in several different units */
  enum shell_units {not_known=0, comoving_distance=1, redshift=2, expansion_factor=3};

  FILE *fd = fopen(radius_file, "r");
  if(!fd)error("Failed to open lightcone radius file %s", radius_file);

  /* Count number of non-zero length lines */
  size_t len = 0;
  char *line = NULL;
  int nr_lines = 0;
  while (getline(&line, &len, fd) != -1 && strlen(line) > 0) nr_lines+=1;
  rewind(fd);

  /* Allocate output array */
  struct lightcone_shell *shell = malloc(sizeof(struct lightcone_shell)*(nr_lines-1));

  /* Check header */
  enum shell_units units = not_known;
  if(getline(&line, &len, fd) != -1) {
    if (strcmp(line, "# Minimum comoving distance, Maximum comoving distance\n") == 0) {
      units = comoving_distance;
    } else if (strcmp(line, "# Minimum redshift, Maximum redshift\n") == 0) {
      units = redshift;
    } else if (strcmp(line, "# Maximum expansion factor, Minimum expansion factor\n") == 0) {
      units = expansion_factor;
    } else {
      error("Unrecognized header in radius file");
    }
  } else {
    error("Unable to read header in radius file");
  }

  /* Read lines */
  for(int i=0; i<nr_lines-1; i+=1) {
    if(fscanf(fd, "%le, %le\n", &shell[i].rmin, &shell[i].rmax) != 2)
      error("Failed to read line from radius file");
  }
  fclose(fd);
  *nr_shells = nr_lines-1;
  const int nr = *nr_shells;
  free(line);

  /* Convert units */
  switch(units) {
  case comoving_distance:
    /* Input is already comoving distance */
    break;
  case redshift:
    /* Convert redshift to comoving distance */
    for(int i=0; i<nr; i+=1) {
      const double a_at_rmin = 1.0/(1.0+shell[i].rmin);
      shell[i].rmin = cosmology_get_comoving_distance(cosmo, a_at_rmin);
      const double a_at_rmax = 1.0/(1.0+shell[i].rmax);
      shell[i].rmax = cosmology_get_comoving_distance(cosmo, a_at_rmax);
    }
    break;
  case expansion_factor:
    /* Convert expansion factor to comoving distance */
    for(int i=0; i<nr; i+=1) {
      shell[i].rmin = cosmology_get_comoving_distance(cosmo, shell[i].rmin);
      shell[i].rmax = cosmology_get_comoving_distance(cosmo, shell[i].rmax);
    }
    break;
  default:
    error("unknown unit type");
  }

  /* Do some sanity checks on the radii */
  /* All values should be monotonically increasing */
  for(int i=1; i<nr; i+=1) {
    if(shell[i].rmin <= shell[i-1].rmin)error("Minimum radii should be monotonically increasing");
    if(shell[i].rmax <= shell[i-1].rmax)error("Maximum radii should be monotonically increasing");
  }

  /* Maximum radius should be greater than minimum */
  for(int i=0; i<nr; i+=1)
    if(shell[i].rmin >= shell[i].rmax)error("Maximum radius should be greater than minimum");

  /* Shells should not overlap */
  for(int i=1; i<nr; i+=1)
    if(shell[i].rmin < shell[i-1].rmax)error("Shells should not overlap");

  /* Return pointer to array */
  *shell_out = shell;
}




struct lightcone_shell *lightcone_shell_array_init(const struct cosmology *cosmo,
                                                   const char *radius_file, int nr_maps,
                                                   struct lightcone_map_type *map_type,
                                                   int nside, size_t total_nr_pix,
                                                   struct lightcone_particle_type *part_type,
                                                   size_t elements_per_block,
                                                   int *nr_shells_out) {
  /* Read in the shell radii */
  int nr_shells;
  struct lightcone_shell *shell;
  if(engine_rank==0)
    read_shell_radii(cosmo, radius_file, &nr_shells, &shell);
#ifdef WITH_MPI
  MPI_Bcast(&nr_shells, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(engine_rank!=0)shell = malloc(sizeof(struct lightcone_shell)*nr_shells);
  MPI_Bcast(shell, sizeof(struct lightcone_shell)*nr_shells, MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

  /* Compute expansion factor at shell edges */
  for(int shell_nr=0; shell_nr<nr_shells; shell_nr+=1) {
    /* Inner edge of the shell */
    shell[shell_nr].amax =
      cosmology_scale_factor_at_comoving_distance(cosmo, shell[shell_nr].rmin);
    /* Outer edge of the shell */
    shell[shell_nr].amin =
      cosmology_scale_factor_at_comoving_distance(cosmo, shell[shell_nr].rmax);
  }

  /* Set initial state of the lightcone shells */
  for(int shell_nr=0;shell_nr<nr_shells; shell_nr+=1)
    shell[shell_nr].state = shell_uninitialized;

  /* Allocate lightcone_map structs for each shell */
  for(int shell_nr=0; shell_nr<nr_shells; shell_nr+=1) {
    shell[shell_nr].nr_maps = nr_maps;
    shell[shell_nr].map = malloc(nr_maps*sizeof(struct lightcone_map));
  }

  /* Initialize lightcone_maps for each shell */
  for(int shell_nr=0;shell_nr<nr_shells; shell_nr+=1) {
    for(int map_nr=0; map_nr<nr_maps; map_nr+=1) {
      lightcone_map_init(&shell[shell_nr].map[map_nr], total_nr_pix,
                         nside, shell[shell_nr].rmin, shell[shell_nr].rmax,
                         map_type[map_nr].units);
    }
  }

  /* Initialize data buffers for map updates - one per particle type per shell */
  for(int shell_nr=0;shell_nr<nr_shells; shell_nr+=1) {
    for(int ptype=0; ptype<swift_type_count; ptype+=1) {
      particle_buffer_init(&shell[shell_nr].buffer[ptype],
                           part_type[ptype].buffer_element_size,
                           elements_per_block, "lightcone_map_updates");
    }
  }

  /* Return the array of shells */
  *nr_shells_out = nr_shells;
  return shell;
}


void lightcone_shell_array_free(struct lightcone_shell *shell, int nr_shells) {

  /* Free the lightcone healpix maps for each shell */
  for(int shell_nr=0;shell_nr<nr_shells; shell_nr+=1) {
    const int nr_maps = shell[shell_nr].nr_maps;
    for(int map_nr=0; map_nr<nr_maps; map_nr+=1) {
      lightcone_map_clean(&shell[shell_nr].map[map_nr]);
    }
  }

  /* Free the arrays of lightcone_map structs for each shell */
  for(int shell_nr=0; shell_nr<nr_shells; shell_nr+=1) {
    free(shell[shell_nr].map);
  }
  
  /* Free the buffers associated with each shell */  
  for(int shell_nr=0;shell_nr<nr_shells; shell_nr+=1) {
    for(int ptype=0; ptype<swift_type_count; ptype+=1) {
      particle_buffer_free(&shell[shell_nr].buffer[ptype]);
    }
  }

  /* Free the array of shells */
  free(shell);

}

void lightcone_shell_array_dump(const struct lightcone_shell *shell, int nr_shells, FILE *stream) {
  
  /* Dump the array of shell structs  */
  restart_write_blocks((void *) shell, sizeof(struct lightcone_shell),
                       nr_shells, stream, "lightcone_shells", "lightcone_shells");

  /* Dump the lightcone maps associated with each shell */
  for(int shell_nr=0;shell_nr<nr_shells; shell_nr+=1) {
    const int nr_maps = shell[shell_nr].nr_maps;
    for(int map_nr=0; map_nr<nr_maps; map_nr+=1) {
      lightcone_map_struct_dump(&shell[shell_nr].map[map_nr], stream);
    }
  }

}


struct lightcone_shell *lightcone_shell_array_restore(FILE *stream, int nr_shells,
                                                      struct lightcone_particle_type *part_type,
                                                      size_t elements_per_block) {

  /* Restore the array of lightcone_shell structs */
  struct lightcone_shell *shell = malloc(sizeof(struct lightcone_shell)*nr_shells);
  restart_read_blocks((void *) shell, sizeof(struct lightcone_shell),
                      nr_shells, stream, NULL, "lightcone_shells");

  /* Restore the lightcone maps associated with each shell */
  for(int shell_nr=0;shell_nr<nr_shells; shell_nr+=1) {
    const int nr_maps = shell[shell_nr].nr_maps;
    shell[shell_nr].map = malloc(sizeof(struct lightcone_map)*nr_maps);
    for(int map_nr=0; map_nr<nr_maps; map_nr+=1) {
      lightcone_map_struct_restore(&shell[shell_nr].map[map_nr], stream);
    }
  }

  /* Initialise the map update buffers */
  for(int shell_nr=0;shell_nr<nr_shells; shell_nr+=1) {
    for(int ptype=0; ptype<swift_type_count; ptype+=1) {
      particle_buffer_init(&shell[shell_nr].buffer[ptype],
                           part_type[ptype].buffer_element_size,
                           elements_per_block, "lightcone_map_updates");
    }
  }

  return shell;
}
