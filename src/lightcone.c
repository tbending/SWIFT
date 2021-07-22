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
#include <math.h>
#include <string.h>
#include <hdf5.h>

/* This object's header. */
#include "lightcone.h"

/* Local headers */
#include "common_io.h"
#include "cosmology.h"
#include "engine.h"
#include "error.h"
#include "extra_io.h"
#include "hydro.h"
#include "lightcone_particle_io.h"
#include "lightcone_replications.h"
#include "lock.h"
#include "parser.h"
#include "part_type.h"
#include "particle_buffer.h"
#include "periodic.h"
#include "restart.h"
#include "space.h"
#include "timeline.h"

/* Whether to dump the replication list */
//#define DUMP_REPLICATIONS
#ifdef DUMP_REPLICATIONS
static int output_nr = 0;
#endif

/* MPI rank for diagnostic messages */
extern int engine_rank;


/**
 * @brief Identify which healpix map types we're making
 *
 * @param props the #lightcone_props structure
 *
 * For each requested map type find the update functions by matching names.
 * Map types are defined in lightcone_map_types.h and there may be extra
 * types defined in extra_io.h.
 *
 * This function assumes that props->map_type is already allocated and
 * props->map_type[:].name has been set to the list of map names from the
 * .yml file. It sets the update_map, ptype_contributes and units fields
 * in the props->map_type array.
 *
 */
static void lightcone_identify_map_types(struct lightcone_props *props) {
  
  /* Loop over requested map types */
  for(int map_nr=0; map_nr<props->nr_maps; map_nr+=1) {

    /* Use null function pointer to indicate not found yet */
    props->map_type[map_nr].update_map = NULL;

    /* First we search the default set of map types in lightcone_map_types.c */
    int type_nr = 0;
    while(lightcone_map_types[type_nr].update_map) {
      if(strcmp(lightcone_map_types[type_nr].name, props->map_type[map_nr].name)==0) {
        props->map_type[map_nr].update_map = lightcone_map_types[type_nr].update_map;
        props->map_type[map_nr].ptype_contributes = lightcone_map_types[type_nr].ptype_contributes;
        props->map_type[map_nr].units = lightcone_map_types[type_nr].units;
        if(engine_rank==0)message("lightcone %d: lightcone map %d is of type %s", 
                                  props->index, map_nr, lightcone_map_types[type_nr].name);
      }
      type_nr += 1;
    }

    /* Then try any additional map types from extra_io.h */
    type_nr = 0;
    while(extra_lightcone_map_types[type_nr].update_map) {
      if(strcmp(extra_lightcone_map_types[type_nr].name, props->map_type[map_nr].name)==0) {
        props->map_type[map_nr].update_map = extra_lightcone_map_types[type_nr].update_map;
        props->map_type[map_nr].ptype_contributes = extra_lightcone_map_types[type_nr].ptype_contributes;
        props->map_type[map_nr].units = extra_lightcone_map_types[type_nr].units;
        if(engine_rank==0)message("lightcone %d: lightcone map %d is of type %s", 
                                  props->index, map_nr, extra_lightcone_map_types[type_nr].name);
      }
      type_nr += 1;
    }

    if(!props->map_type[map_nr].update_map)error("Unable to locate lightcone map type %s",
                                                 props->map_type[map_nr].name);
  }
}


/**
 * @brief Allocate particle I/O buffers for a lightcone
 *
 * @param props the #lightcone_props structure
 */
static void lightcone_allocate_buffers(struct lightcone_props *props) {

  /* Initialize particle output buffers */
  const size_t elements_per_block = (size_t) props->buffer_chunk_size;

  if(props->use_type[swift_type_gas]) {
    particle_buffer_init(&props->buffer[swift_type_gas],
                         sizeof(struct lightcone_gas_data),
                         elements_per_block, "lightcone_gas");
  }
  
  if(props->use_type[swift_type_dark_matter]) {
    particle_buffer_init(&props->buffer[swift_type_dark_matter],
                         sizeof(struct lightcone_dark_matter_data),
                         elements_per_block, "lightcone_dm");
  }

  if(props->use_type[swift_type_dark_matter_background]) {  
    particle_buffer_init(&props->buffer[swift_type_dark_matter_background],
                         sizeof(struct lightcone_dark_matter_data),
                         elements_per_block, "lightcone_dm_bg");
  }

  if(props->use_type[swift_type_stars]) {  
    particle_buffer_init(&props->buffer[swift_type_stars],
                         sizeof(struct lightcone_stars_data),
                         elements_per_block, "lightcone_stars");
  }

  if(props->use_type[swift_type_black_hole]) {
  particle_buffer_init(&props->buffer[swift_type_black_hole],
                       sizeof(struct lightcone_black_hole_data),
                       elements_per_block, "lightcone_bh");
  }

  if(props->use_type[swift_type_neutrino]) {
    particle_buffer_init(&props->buffer[swift_type_neutrino],
                         sizeof(struct lightcone_neutrino_data),
                         elements_per_block, "lightcone_neutrino");
  }  
}


/**
 * @brief Dump lightcone_props struct to the output stream.
 *
 * @param props the #lightcone_props structure
 * @param stream The stream to write to.
 */
void lightcone_struct_dump(const struct lightcone_props *props, FILE *stream) {

  /* Don't dump the replication list - will regenerate it as needed */
  struct lightcone_props tmp = *props;
  tmp.replication_list.nrep = 0;
  tmp.replication_list.replication = NULL;
  tmp.have_replication_list = 0;

  /* Don't write out particle buffers - must flush before dumping restart. */
  memset(tmp.buffer, 0, sizeof(struct particle_buffer)*swift_type_count);

  /* Don't write array pointers */
  tmp.shell = NULL;
  tmp.map_type = NULL;
  for(int ptype=0; ptype<swift_type_count; ptype+=1)
    tmp.part_type[ptype].map_index = NULL;

  /* Dump the lightcone struct */
  restart_write_blocks((void *) &tmp, sizeof(struct lightcone_props), 1, stream,
                       "lightcone_props", "lightcone_props");

  /* Dump the array of map types */
  restart_write_blocks((void *) props->map_type, sizeof(struct lightcone_map_type),
                       props->nr_maps, stream, "lightcone_props", "lightcone_props");

  /* Dump the array of shells */
  lightcone_shell_array_dump(props->shell, props->nr_shells, stream);

  /* For each particle type we have an array of lightcone map indexes to update. Dump these. */
  for(int ptype=0; ptype<swift_type_count; ptype+=1) {
    const struct lightcone_particle_type *this_type = &(props->part_type[ptype]);
    restart_write_blocks((void *) this_type->map_index, sizeof(int),
                         this_type->nr_maps, stream, "lightcone_props", "lightcone_props");
  }
}


/**
 * @brief Restore lightcone_props struct from the input stream.
 *
 * @param props the #lightcone_props structure
 * @param stream The stream to read from.
 */
void lightcone_struct_restore(struct lightcone_props *props, FILE *stream) {

  /* Restore lightcone struct */
  restart_read_blocks((void *)props, sizeof(struct lightcone_props), 1, stream,
                      NULL, "lightcone_props");

  /* Read in the map types */
  props->map_type = malloc(sizeof(struct lightcone_map_type)*props->nr_maps);
  restart_read_blocks((void *)props->map_type, sizeof(struct lightcone_map_type),
                      props->nr_maps, stream, NULL, "lightcone_props");

  /* Read in the shells */
  props->shell = lightcone_shell_array_restore(stream, props->nr_shells, props->part_type,
                                               props->buffer_chunk_size);

  /* For each particle type we have an array of lightcone map indexes to update. Restore these. */
  for(int ptype=0; ptype<swift_type_count; ptype+=1) {
    struct lightcone_particle_type *this_type = &(props->part_type[ptype]);
    this_type->map_index = malloc(sizeof(int)*this_type->nr_maps);
    restart_read_blocks((void *)this_type->map_index, sizeof(int),
                        this_type->nr_maps, stream, NULL, "lightcone_props");
  }

  /* Restore pointers to functions for updating healpix maps */
  lightcone_identify_map_types(props);

  /* Re-initialise C++ smoothing code */
  props->smoothing_info = healpix_smoothing_init(props->nside, kernel_gamma);

  /* Re-allocate particle data buffers */
  lightcone_allocate_buffers(props);

  /* Define output quantities */
  lightcone_io_make_output_fields();
}


static char *yaml_name(char *buf, const char *str1, const char *str2) {
  int len = snprintf(buf, PARSER_MAX_LINE_SIZE, "%s:%s", str1, str2);
  if((len < 0) || (len >= PARSER_MAX_LINE_SIZE))
    error("Failed to generate parameter name");
  return buf;
}


/**
 * @brief Initialise the properties of the lightcone code.
 *
 * @param props the #lightcone_props structure to fill.
 * @param s the #space structure.
 * @param cosmo the #cosmology structure.
 * @param params the parameter file parser.
 * @param verbose the verbosity flag
 */
void lightcone_init(struct lightcone_props *props,
                    const char *name, int index,
                    const struct space *s,
                    const struct cosmology *cosmo,
                    struct swift_params *params,
                    const int verbose) {
  
  /* Macro to generate parameter names given section name */
  char buf[PARSER_MAX_LINE_SIZE];
#define YML_NAME(x) yaml_name(buf, name, x)

  /* Store index of this lightcone in the .yml file */
  props->index = index;

  /* Verbose lightcone output - use passed in value of --verbose flag */
  props->verbose = verbose;

  /* Define output quantities */
  lightcone_io_make_output_fields();

  /* Which particle types we should write out particle data for */
  for(int i=0; i<swift_type_count; i+=1)
    props->use_type[i] = 0;
  props->use_type[swift_type_gas] = parser_get_param_int(params, YML_NAME("use_gas"));
  props->use_type[swift_type_dark_matter] = parser_get_param_int(params, YML_NAME("use_dm"));
  props->use_type[swift_type_dark_matter_background] = parser_get_param_int(params, YML_NAME("use_dm_background"));
  props->use_type[swift_type_stars] = parser_get_param_int(params, YML_NAME("use_stars"));
  props->use_type[swift_type_black_hole] = parser_get_param_int(params, YML_NAME("use_black_hole"));
  props->use_type[swift_type_neutrino] = parser_get_param_int(params, YML_NAME("use_neutrino"));

  /* Base name for output files */
  parser_get_param_string(params, YML_NAME("basename"), props->basename);

  /* Redshift range for particle output */
  props->z_min_for_particles = parser_get_param_double(params, YML_NAME("z_min_for_particles"));
  props->z_max_for_particles = parser_get_param_double(params, YML_NAME("z_max_for_particles"));

  /* Corresponding range in comoving distance squared */
  const double a_min_for_particles = 1.0/(1.0+props->z_max_for_particles);
  props->r2_max_for_particles = pow(cosmology_get_comoving_distance(cosmo, a_min_for_particles), 2.0);
  const double a_max_for_particles = 1.0/(1.0+props->z_min_for_particles);
  props->r2_min_for_particles = pow(cosmology_get_comoving_distance(cosmo, a_max_for_particles), 2.0);

  /* Coordinates of the observer in the simulation box */
  parser_get_param_double_array(params, YML_NAME("observer_position"), 3,
                                props->observer_position);

  /* Write particles to disk if this many or more are in the buffer */
  props->max_particles_buffered = parser_get_opt_param_int(params, YML_NAME("max_particles_buffered"), 100000);

  /* Chunk size for particles buffered in memory  */
  props->buffer_chunk_size = parser_get_opt_param_int(params, YML_NAME("buffer_chunk_size"), 20000);

  /* Chunk size for particles buffered in memory  */
  props->hdf5_chunk_size = parser_get_opt_param_int(params, YML_NAME("hdf5_chunk_size"), 16384);

  /* Get the size of the simulation box */
  props->boxsize = s->dim[0];
  if(s->dim[1] != s->dim[0] || s->dim[2] != s->dim[0])
    error("Lightcones require a cubic simulation box.");

  /* Get top level cell size */
  props->cell_width = s->width[0];
  if(s->width[1] != s->width[0] || s->width[2] != s->width[0])
    error("Lightcones require cubic top level cells.");

  /* Initially have no replication list */
  props->have_replication_list = 0;
  props->ti_old = 0;
  props->ti_current = 0;

  /* Initialize various counters */
  for(int i=0; i<swift_type_count; i+=1) {
    props->tot_num_particles_written[i] = 0;
    props->num_particles_written_to_file[i] = 0;
  }
  props->current_file = -1;
  
  /* Always start a new file initially */
  props->start_new_file = 1;

  /* 
     Healpix map parameters for this lightcone
  */

  /* Healpix nside parameter */
  props->nside = parser_get_param_double(params, YML_NAME("nside"));

  /* Initialise C++ smoothing code */
  props->smoothing_info = healpix_smoothing_init(props->nside, kernel_gamma);

  /* Update lightcone pixel data if more than this number of updates are buffered */
  props->max_updates_buffered = parser_get_opt_param_int(params, YML_NAME("max_updates_buffered"), 1000000);
  
  /* Name of the file with radii of spherical shells */
  parser_get_param_string(params, YML_NAME("radius_file"), props->radius_file);
  
  /* Whether we smooth the maps */
  props->smooth = parser_get_opt_param_int(params, YML_NAME("smooth"), 1);

  /* Names of the healpix maps to make for this lightcone */
  char **map_names;
  parser_get_param_string_array(params, YML_NAME("map_names"), &props->nr_maps, &map_names);
  props->map_type = malloc(props->nr_maps*sizeof(struct lightcone_map_type));
  for(int i=0; i<props->nr_maps; i+=1) {
    int len = snprintf(props->map_type[i].name, PARSER_MAX_LINE_SIZE, "%s", map_names[i]);
    if(len >= PARSER_MAX_LINE_SIZE || len < 0) error("Map type name truncated or encoding error");
  }
  parser_free_param_string_array(props->nr_maps, map_names);

  /* For each requested map type find the update function by matching names */
  lightcone_identify_map_types(props);

  /* For each particle type, determine which healpix maps will be updated */
  const int nr_maps = props->nr_maps;
  for(int ptype=0; ptype<swift_type_count; ptype+=1) {
    
    struct lightcone_particle_type *this_type = &(props->part_type[ptype]);

    /* Count maps updated by this particle type */
    this_type->nr_maps = 0;
    for(int map_nr=0; map_nr<nr_maps; map_nr+=1) {
      if(props->map_type[map_nr].ptype_contributes(ptype))
        this_type->nr_maps += 1;
    }
    
    /* Store indexes of maps to update for this particle type */
    this_type->map_index = malloc(sizeof(int)*this_type->nr_maps);
    this_type->nr_maps = 0;
    for(int map_nr=0; map_nr<nr_maps; map_nr+=1) {
      if(props->map_type[map_nr].ptype_contributes(ptype)) {
        this_type->map_index[this_type->nr_maps] = map_nr;
        this_type->nr_maps += 1;
      }
    }

    /* Determine how much data we need to store per particle of this type.
       We need theta and phi angular coordinates, angular size of the particle,
       and the values to be added to the healpix maps */
    this_type->buffer_element_size = (3+this_type->nr_maps) * sizeof(double);
  }

  /* Set up the array of lightcone shells for this lightcone */
  size_t total_nr_pix = healpix_smoothing_get_npix(props->smoothing_info);
  props->shell = lightcone_shell_array_init(cosmo, props->radius_file, props->nr_maps,
                                            props->map_type, props->nside, total_nr_pix,
                                            props->part_type, props->buffer_chunk_size,
                                            &props->nr_shells);

  /* Compute area of a healpix pixel */
  props->pixel_area_steradians = 4*M_PI/total_nr_pix;
  
  /* Report shell radii */
  const int nr_shells = props->nr_shells;
  if(engine_rank==0) {
    for(int shell_nr=0; shell_nr<nr_shells; shell_nr+=1) {
      message("lightcone %d: shell %d has inner radius %e and outer radius %e", 
              index, shell_nr, props->shell[shell_nr].rmin, props->shell[shell_nr].rmax);
    }
  }
  if(engine_rank==0)message("lightcone %d: there are %d lightcone shells and %d maps per shell",
                            index, nr_shells, nr_maps);

  /* Determine full redshift range to search for lightcone crossings.
     Find range in expansion factor for particle output. */
  double a_min = 1.0/(1.0+props->z_max_for_particles);
  double a_max = 1.0/(1.0+props->z_min_for_particles);
  /* Then extend the range to include all healpix map shells */
  for(int shell_nr=0;shell_nr<nr_shells; shell_nr+=1) {
    const double shell_a_min = props->shell[shell_nr].amin;
    const double shell_a_max = props->shell[shell_nr].amax;
    if(shell_a_min < a_min)a_min = shell_a_min;
    if(shell_a_max > a_max)a_max = shell_a_max;
  }
  props->a_min = a_min;
  props->a_max = a_max;
  if(engine_rank==0)message("lightcone %d: range in expansion factor to search lightcone: %e to %e",
                            index, a_min, a_max);

  /* Store the corresponding comoving distance squared */
  props->r2_max = pow(cosmology_get_comoving_distance(cosmo, a_min), 2.0);
  props->r2_min = pow(cosmology_get_comoving_distance(cosmo, a_max), 2.0);

  /* Allocate lightcone output buffers */
  lightcone_allocate_buffers(props);

  /* Estimate number of particles which will be output.
     
     Assumptions:
     - flat cosmology (haven't implemented comoving volume calculation for non-flat)
     - uniform box
  */
  const long long nr_gparts = s->nr_gparts;
  long long total_nr_gparts;
#ifdef WITH_MPI
  MPI_Reduce(&nr_gparts, &total_nr_gparts, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
#else
  total_nr_gparts = nr_gparts;
#endif
  if(engine_rank==0) {
    const double lightcone_rmax = cosmology_get_comoving_distance(cosmo, a_min_for_particles);
    const double lightcone_rmin = cosmology_get_comoving_distance(cosmo, a_max_for_particles);
    const double volume = 4./3.*M_PI*(pow(lightcone_rmax, 3.)-pow(lightcone_rmin, 3.));
    const long long est_nr_output = total_nr_gparts / pow(props->boxsize, 3.0) * volume;
    message("lightcone %d: comoving distance to max. particle redshift: %e", index, lightcone_rmax);
    message("lightcone %d: gparts in lightcone (if uniform box+flat cosmology): %lld", index, est_nr_output);
  }

}


static void particle_file_name(char *buf, int len, char *basename,
                               int current_file, int comm_rank) {
  
  int ret = snprintf(buf, len, "%s_%04d.%d.hdf5", basename, current_file, comm_rank);
  if((ret < 0) || (ret >= len))error("Lightcone particle file name truncation or output error");
}


/**
 * @brief Flush any buffers which exceed the specified size.
 *
 * @param props the #lightcone_props structure.
 * @param flush_all flag to force flush of all buffers
 * @param end_file if true, subsequent calls write to a new file
 *
 */
void lightcone_flush_particle_buffers(struct lightcone_props *props,
                                      const struct unit_system *internal_units,
                                      const struct unit_system *snapshot_units,
                                      int flush_all, int end_file) {

  ticks tic = getticks();
  
  /* Will flush any buffers with more particles than this */
  size_t max_to_buffer = (size_t) props->max_particles_buffered;
  if(flush_all)max_to_buffer = 0;

  /* Count how many types have data to write out */
  int types_to_flush = 0;
  for(int ptype=0; ptype<swift_type_count; ptype+=1) {
    if(props->use_type[ptype]) {
      const size_t num_to_write = particle_buffer_num_elements(&props->buffer[ptype]);
      if(num_to_write >= max_to_buffer && num_to_write > 0)types_to_flush += 1;
    }
  }
  
  /* Check if there's anything to do */
  if(types_to_flush>0) {
    
    /* We have data to flush, so open or create the output file */
    hid_t file_id;
    char fname[FILENAME_BUFFER_SIZE];
    if(props->start_new_file) {

      /* Get the name of the next file */
      props->current_file += 1;
      particle_file_name(fname, FILENAME_BUFFER_SIZE, props->basename,
                         props->current_file, engine_rank);

      /* Create the file */
      file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      if(file_id < 0)error("Unable to create new lightcone file: %s", fname);

      /* We have now written no particles to the current file */
      for(int ptype=0; ptype<swift_type_count; ptype+=1)
        props->num_particles_written_to_file[ptype] = 0;    

      /* Write the system of Units used in the spashot */
      io_write_unit_system(file_id, snapshot_units, "Units");

      /* Write the system of Units used internally */
      io_write_unit_system(file_id, internal_units, "InternalCodeUnits");

      /* Write the observer position and redshift limits */
      hid_t group_id = H5Gcreate2(file_id, "Lightcone", H5P_DEFAULT,
                                  H5P_DEFAULT, H5P_DEFAULT);
      io_write_attribute(group_id, "observer_position", DOUBLE,
                         props->observer_position, 3);
      io_write_attribute_d(group_id, "minimum_redshift", props->z_min_for_particles);
      io_write_attribute_d(group_id, "maximum_redshift", props->z_max_for_particles);

      /* Record number of MPI ranks so we know how many files there are */
      int comm_rank = 0;
      int comm_size = 1;
#ifdef WITH_MPI
      MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
      MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif
      io_write_attribute_i(group_id, "mpi_rank", comm_rank);
      io_write_attribute_i(group_id, "nr_mpi_ranks", comm_size);
      io_write_attribute_i(group_id, "file_index", props->current_file);

      H5Gclose(group_id);

      /* We no longer need to create a new file */
      props->start_new_file = 0;

    } else {

      /* Re-open an existing file */
      particle_file_name(fname, FILENAME_BUFFER_SIZE, props->basename,
                         props->current_file, engine_rank);
      file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);
      if(file_id < 0)error("Unable to open current lightcone file: %s", fname);

    }

    /* Loop over particle types */
    for(int ptype=0; ptype<swift_type_count; ptype+=1) {
      if(props->use_type[ptype]) {
        const size_t num_to_write = particle_buffer_num_elements(&props->buffer[ptype]);
        if(num_to_write >= max_to_buffer && num_to_write > 0) {
          if(props->verbose)message("lightcone %d: dumping %d particles of type %s",
                                    props->index, (int) num_to_write, part_type_names[ptype]);
          lightcone_write_particles(props, internal_units, snapshot_units, ptype, file_id);
          particle_buffer_empty(&props->buffer[ptype]);
          props->num_particles_written_to_file[ptype] += num_to_write;          
        }
      }
    }

    /* We're done updating the output file */
    H5Fclose(file_id);
  }

  /* If we need to start a new file next time, record this */
  if(end_file)props->start_new_file = 1;

  if (props->verbose && engine_rank == 0 && types_to_flush > 0)
    message("lightcone %d: Flushing particle buffers took %.3f %s.",
            props->index, clocks_from_ticks(getticks() - tic), clocks_getunit());

}


/**
 * @brief Flush lightcone map update buffers for all shells
 *
 * @param props the #lightcone_props structure.
 *
 */
void lightcone_flush_map_updates(struct lightcone_props *props,
                                 struct threadpool *tp) {    

  ticks tic = getticks();

  /* Report how much memory we're using before flushing buffers */
  if(props->verbose)lightcone_report_memory_use(props);

  /* Apply updates to all current shells */
  for(int shell_nr=0; shell_nr<props->nr_shells; shell_nr+=1) {
    if(props->shell[shell_nr].state == shell_current)
      lightcone_shell_flush_map_updates(&props->shell[shell_nr], tp,
                                        props->part_type, props->smoothing_info);
  }

  /* Report runtime */
  if (props->verbose && engine_rank==0)
    message("lightcone %d: Applying lightcone map updates took %.3f %s.",
            props->index, clocks_from_ticks(getticks() - tic), clocks_getunit());
}


/**
 * @brief Write and deallocate any completed lightcone shells
 *
 * @param props the #lightcone_props structure.
 * @param c the #cosmology structure
 * @param dump_all flag to indicate that all shells should be dumped
 * @param need_flush whether there might be buffered updates to apply
 *
 */
void lightcone_dump_completed_shells(struct lightcone_props *props,
                                     struct threadpool *tp,
                                     const struct cosmology *c,
                                     const struct unit_system *internal_units,
                                     const struct unit_system *snapshot_units,
                                     const int dump_all,
                                     const int need_flush) {
#ifdef HAVE_HDF5

  ticks tic = getticks();

  /* Get number of shells and maps per shell */
  const int nr_shells = props->nr_shells;
  const int nr_maps = props->nr_maps;

  /* Compute expansion factor corresponding to time props->ti_old,
     which is the earliest time any particle might have been drifted
     from on this step. Here we assume that no particle remains to
     be drifted from any time earlier than this so that any shell
     whose redshift range is entirely before ti_old can be now be
     written out and deallocated. */
  const double a_complete = c->a_begin * exp(props->ti_old * c->time_base);

  int num_shells_written = 0;
  for(int shell_nr=0; shell_nr<nr_shells; shell_nr+=1) {

    /* Will write out this shell if it has been updated but not written
       out yet and either we advanced past its redshift range or we're
       dumping all remaining shells at the end of the simulation */
    if(props->shell[shell_nr].state==shell_current) {
      if(props->shell[shell_nr].amax < a_complete || dump_all) {

        if(props->verbose && engine_rank==0)
          message("lightcone %d: writing out completed shell %d at a=%f",
                  props->index, shell_nr, c->a);

        if(num_shells_written==0) {
          /* Report how much memory we're using before flushing buffers */
          lightcone_report_memory_use(props);
        }

        num_shells_written += 1;

        /* Apply any buffered updates for this shell, if we didn't already */
        if(need_flush) lightcone_shell_flush_map_updates(&props->shell[shell_nr], tp,
                                                         props->part_type, props->smoothing_info);

        /* Get the name of the file to write */
        char fname[FILENAME_BUFFER_SIZE];
        int len = snprintf(fname, FILENAME_BUFFER_SIZE, "%s.shell_%d.hdf5",
                           props->basename, shell_nr);
        if((len < 0) || (len >= FILENAME_BUFFER_SIZE))
          error("Lightcone map output filename truncation or output error");
        
        /* Create the output file for this shell */
        hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
#ifdef WITH_MPI
#ifdef HAVE_PARALLEL_HDF5
        if(H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL) < 0)
          error("Unable to set HDF5 MPI-IO file access mode");
#else
        error("Writing lightcone maps with MPI requires parallel HDF5");
#endif
#endif
        hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
        if(file_id < 0)error("Unable to create file %s", fname);

        /* Write the system of Units used in the spashot */
        io_write_unit_system(file_id, snapshot_units, "Units");

        /* Write the system of Units used internally */
        io_write_unit_system(file_id, internal_units, "InternalCodeUnits");

        /* Write the lightcone maps for this shell */
        for(int map_nr=0; map_nr<nr_maps; map_nr+=1)
          lightcone_map_write(&(props->shell[shell_nr].map[map_nr]), file_id, props->map_type[map_nr].name,
                              internal_units, snapshot_units);

        /* Close the file */
        H5Pclose(fapl_id);
        H5Fclose(file_id);

        /* Free the pixel data associated with this shell */
        for(int map_nr=0; map_nr<nr_maps; map_nr+=1)
          lightcone_map_free_pixels(&(props->shell[shell_nr].map[map_nr]));

        /* Update status of this shell */
        props->shell[shell_nr].state=shell_complete;
      }
    }
  }

  if (props->verbose && engine_rank==0 && num_shells_written > 0)
    message("lightcone %d: Writing completed lightcone shells took %.3f %s.",
            props->index, clocks_from_ticks(getticks() - tic), clocks_getunit());

#else
  error("Need HDF5 to write out lightcone maps");
#endif
}


/**
 * @brief Deallocate lightcone data.
 *
 * @param props the #lightcone_props structure.
 *
 */
void lightcone_clean(struct lightcone_props *props) {

  /* Deallocate particle buffers */
  for(int i=0; i<swift_type_count; i+=1) {
    if(props->use_type[i])
      particle_buffer_free(&props->buffer[i]);
  }

  /* Free replication list, if we have one */
  if(props->have_replication_list)replication_list_clean(&props->replication_list);

  /* Clean lightcone maps and free the structs */
  const int nr_shells = props->nr_shells;
  const int nr_maps = props->nr_maps;
  for(int shell_nr=0;shell_nr<nr_shells; shell_nr+=1) {
    for(int map_nr=0; map_nr<nr_maps; map_nr+=1) {
      lightcone_map_clean(&(props->shell[shell_nr].map[map_nr]));
    }
    free(props->shell[shell_nr].map);
  }

  /* Free buffers associated with each shell */
  for(int shell_nr=0;shell_nr<nr_shells; shell_nr+=1) {
    for(int ptype=0; ptype<swift_type_count; ptype+=1) {
      particle_buffer_free(&(props->shell[shell_nr].buffer[ptype]));
    }
  }
  
  /* Free array of shells */
  free(props->shell);

  /* Free array of lightcone map types */
  free(props->map_type);

  /* Free data associated with particle types */
  for(int ptype=0; ptype<swift_type_count; ptype+=1) {
    struct lightcone_particle_type *this_type = &(props->part_type[ptype]);
    free(this_type->map_index);
  }
  
  /* Tidy up C++ healpix smoothing code */
  healpix_smoothing_clean(props->smoothing_info);
}


/**
 * @brief Determine periodic copies of the simulation box which could
 * contribute to the lightcone.
 *
 *                     \
 *           \          \
 *            |         |
 * Obs      A |    B    | C
 *            |         |
 *           /          /
 *          R1         /
 *                    R0
 *
 * Consider a single particle being drifted. Here R0 is the comoving
 * distance to the time the particle is drifted FROM. R1 is the comoving
 * distance to the time the particle is drifted TO on this step.
 *
 * Particles which are beyond the lightcone surface at the start of
 * their drift (C) cannot cross the lightcone on this step if v < c.
 * Particles between the lightcone surfaces at the start and end of
 * their drift (B) may cross the lightcone (and certainly will if they
 * have zero velocity).
 *
 * Particles just within the lightcone surface at the start of their
 * drift (A) may be able to cross the lightcone due to their velocity so
 * we need to allow a boundary layer on the inside edge of the shell.
 * If we assume v < c, then we can use a layer of thickness R0-R1.
 *
 * Here we compute the earliest and latest times particles may be drifted
 * between, find the corresponding comoving distances R0 and R1, reduce
 * the inner distance by R0-R1, and find all periodic copies of the 
 * simulation box which overlap this spherical shell.
 *
 * Later we use this list to know which periodic copies to check when
 * particles are drifted.
 *
 * This routine also determines which lightcone healpix maps might be
 * updated on this time step and allocates the pixel data if necessary.
 *
 * @param props The #lightcone_props structure
 * @param cosmo The #cosmology structure
 * @param ti_old Beginning of the timestep
 * @param ti_current End of the timestep
 * @param dt_max Maximum time step length
 */
void lightcone_prepare_for_step(struct lightcone_props *props,
                                const struct cosmology *cosmo,
                                const integertime_t ti_old,
                                const integertime_t ti_current,
                                const double dt_max) {
  ticks tic = getticks();

  /* Deallocate the old list, if there is one */
  if(props->have_replication_list)replication_list_clean(&props->replication_list);

  /* Get the size of the simulation box */
  const double boxsize = props->boxsize;

  /* Get a lower limit on earliest time particle may be drifted from */
  float dt = cosmo->time_end - cosmo->time_begin;
  while (dt > dt_max) dt /= 2.f;
  timebin_t bin = get_time_bin(dt*cosmo->time_base_inv);
  integertime_t ti_lim = get_integer_time_begin(ti_old, bin);

  /* Get expansion factor at earliest and latest times particles might be drifted between */
  double a_current = cosmo->a_begin * exp(ti_current * cosmo->time_base);
  double a_old = cosmo->a_begin * exp(ti_lim * cosmo->time_base);
  if(a_old < cosmo->a_begin)a_old = cosmo->a_begin;

  /* Convert redshift range to a distance range */
  double lightcone_rmin = cosmology_get_comoving_distance(cosmo, a_current);
  double lightcone_rmax = cosmology_get_comoving_distance(cosmo, a_old);
  if(lightcone_rmin > lightcone_rmax)
    error("Lightcone has rmin > rmax");

  /* Allow inner boundary layer, assuming all particles have v < c.
     This is to account for particles moving during the time step. */
  double boundary = lightcone_rmax-lightcone_rmin;
  lightcone_rmin -= boundary;
  if(lightcone_rmin < 0)lightcone_rmin = 0;

  if(a_current < props->a_min || a_old > props->a_max) {
    /* Timestep does not overlap the lightcone redshift range */
    replication_list_init_empty(&props->replication_list);
  } else {
    /* Timestep may contribute particles to the lightcone */
    replication_list_init(&props->replication_list, boxsize,
                          props->cell_width,
                          props->observer_position,
                          lightcone_rmin, lightcone_rmax);
  }

  /* Record that we made the list */
  props->have_replication_list = 1;

  /* Store times we used to make the list, for consistency check later */
  props->ti_old = ti_lim;
  props->ti_current = ti_current;

  /* Report the size of the list */
#ifdef DUMP_REPLICATIONS
  if(engine_rank==0) {
    message("lightcone %d: no. of replications to check: %d", props->index, props->replication_list.nrep);
    message("lightcone %d: shell to search inner radius=%e, outer radius=%e", props->index, lightcone_rmin,
            lightcone_rmax);
  }
#endif

  /* Write out the list, if required */
#ifdef DUMP_REPLICATIONS
  if(engine_rank==0) {
    char fname[500];
    sprintf(fname, "replication_list.%d.txt", output_nr);
    FILE *fd_rep = fopen(fname, "w");
    fprintf(fd_rep, "# Observer x, y, z\n");
    fprintf(fd_rep, "%e, %e, %e\n", props->observer_position[0],
            props->observer_position[1], props->observer_position[2]); 
    fprintf(fd_rep, "# Box size, inner radius, outer radius\n");
    fprintf(fd_rep, "%e, %e, %e\n", boxsize, lightcone_rmin-boundary, lightcone_rmax);
    fprintf(fd_rep, "# x, y, z, rmin2, rmax2\n");
    replication_list_write(&props->replication_list, fd_rep);
    fclose(fd_rep);
    output_nr += 1;
  }
#endif

  /* Number of shells and maps per shell */
  const int nr_maps   = props->nr_maps;
  const int nr_shells = props->nr_shells;

  /* Range of shells that might be updated this step */
  int shell_nr_min = nr_shells;
  int shell_nr_max = -1;

  /* Loop over healpix map shells */
  for(int shell_nr=0; shell_nr<nr_shells; shell_nr+=1) {

    const double shell_amin = props->shell[shell_nr].amin;
    const double shell_amax = props->shell[shell_nr].amax;
    const double step_amin = a_old;
    const double step_amax = a_current;

    /* Check if this shell might be updated */
    if(step_amin <= shell_amax && step_amax >= shell_amin) {

      switch(props->shell[shell_nr].state) {
      case shell_uninitialized:
        /* This shell has not been allocated yet, so allocate it */
        if(props->verbose && engine_rank==0)
          message("lightcone %d: allocating pixels for shell %d at a=%f", props->index, shell_nr, cosmo->a);
        for(int map_nr=0; map_nr<nr_maps; map_nr+=1)
          lightcone_map_allocate_pixels(&(props->shell[shell_nr].map[map_nr]), /* zero_pixels = */ 1);   
        props->shell[shell_nr].state = shell_current;
        break;
      case shell_complete:
        /* Shell has already been written out and freed - should never happen */
        error("Lightcone shell has been written out while particles could still contribute");
        break;
      case shell_current:
        /* Already initialized, nothing to do */
        break;
      }

      /* Record range of shells that might be updated this step */
      if (shell_nr < shell_nr_min) shell_nr_min = shell_nr;
      if (shell_nr > shell_nr_max) shell_nr_max = shell_nr;
    }
  }
  props->shell_nr_min = shell_nr_min;
  props->shell_nr_max = shell_nr_max;

  if (props->verbose && engine_rank==0)
    message("lightcone %d: Lightcone timestep preparations took %.3f %s.",
            props->index, clocks_from_ticks(getticks() - tic), clocks_getunit());

}


/**
 * @brief Determine whether lightcone map buffers should be flushed this step.
 *
 * @param props The #lightcone_props structure
 *
 */
int lightcone_trigger_map_update(struct lightcone_props *props) {
  
  size_t total_updates = 0;
  const int nr_shells = props->nr_shells;
  for(int shell_nr=0; shell_nr<nr_shells; shell_nr+=1) {
    if(props->shell[shell_nr].state == shell_current) {
      for(int ptype=0; ptype<swift_type_count; ptype+=1) {
        total_updates += particle_buffer_num_elements(&(props->shell[shell_nr].buffer[ptype]));
      }
    }
  }
  return total_updates >= ((size_t) props->max_updates_buffered);
}


/**
 * @brief Add a particle to the output buffer
 *
 * @param props The #lightcone_props structure
 * @param e The #engine structure
 * @param gp The #gpart to buffer
 * @param a_cross Expansion factor of lightcone crossing
 * @param x_cross Position of the gpart at lightcone crossing
 */
void lightcone_buffer_particle(struct lightcone_props *props,
                               const struct engine *e, const struct gpart *gp,
                               const double a_cross, const double x_cross[3]) {

  /* Handle on the other particle types */
  const struct space *s = e->s;
  const struct part *parts = s->parts;
  const struct xpart *xparts = s->xparts;
  const struct spart *sparts = s->sparts;
  const struct bpart *bparts = s->bparts;

  switch (gp->type) {
  case swift_type_gas: {

    const struct part *p = &parts[-gp->id_or_neg_offset];
    const struct xpart *xp = &xparts[-gp->id_or_neg_offset];
    struct lightcone_gas_data data;
    if(lightcone_store_gas(e, gp, p, xp, a_cross, x_cross, &data))
      particle_buffer_append(props->buffer+swift_type_gas, &data);
 
  } break;

  case swift_type_stars: {
    
    const struct spart *sp = &sparts[-gp->id_or_neg_offset];
    struct lightcone_stars_data data;
    if(lightcone_store_stars(e, gp, sp, a_cross, x_cross, &data))
      particle_buffer_append(props->buffer+swift_type_stars, &data);

  } break;

  case swift_type_black_hole: {
      
    const struct bpart *bp = &bparts[-gp->id_or_neg_offset];
    struct lightcone_black_hole_data data;
    if(lightcone_store_black_hole(e, gp, bp, a_cross, x_cross, &data))
      particle_buffer_append(props->buffer+swift_type_black_hole, &data);

  } break;

  case swift_type_dark_matter: {

    struct lightcone_dark_matter_data data;
    if(lightcone_store_dark_matter(e, gp, a_cross, x_cross, &data))
      particle_buffer_append(props->buffer+swift_type_dark_matter, &data);
    
  } break;

  case swift_type_dark_matter_background: {

    /* Assumed to have same properties as DM particles */
    struct lightcone_dark_matter_data data;
    if(lightcone_store_dark_matter(e, gp, a_cross, x_cross, &data))
      particle_buffer_append(props->buffer+swift_type_dark_matter_background, &data);
    
  } break;

  case swift_type_neutrino: {

    struct lightcone_neutrino_data data;
    if(lightcone_store_neutrino(e, gp, a_cross, x_cross, &data))
      particle_buffer_append(props->buffer+swift_type_neutrino, &data);

  } break;

  default:
    error("Particle type not supported in lightcones");
  }
  
}


static double angular_smoothing_scale(const double *pos, const double hsml) {
  
  /* Compute distance to particle */
  double dist = 0;
  for(int i=0; i<3; i+=1)
    dist += pos[i]*pos[i];
  dist = sqrt(dist);
  
  /* Avoid trig call for small angles (accurate to about 0.3%) */
  if(dist > 10.0*hsml)
    return hsml/dist;
  else
    return atan(hsml/dist);
}


/**
 * @brief Buffer a particle's contribution to the healpix map(s)
 *
 * @param props The #lightcone_props structure
 * @param e The #engine structure
 * @param gp The #gpart to buffer
 * @param a_cross Expansion factor of lightcone crossing
 * @param x_cross Position of the gpart at lightcone crossing
 *
 */
void lightcone_buffer_map_update(struct lightcone_props *props,
                                 const struct engine *e, const struct gpart *gp,
                                 const double a_cross, const double x_cross[3]) {

  /* Find information on healpix maps this particle type contributes to */
  const struct lightcone_particle_type *part_type_info = &(props->part_type[gp->type]);

  /* If this particle type contributes to no healpix maps, do nothing */
  if(part_type_info->nr_maps==0)return;

  /* Get angular coordinates of the particle */
  double theta, phi;
  healpix_smoothing_vec2ang(props->smoothing_info, x_cross, &theta, &phi);

  /* Get angular size of the particle */
  double radius;
  if((gp->type == swift_type_gas) && props->smooth) {
    const struct part *parts = e->s->parts;
    const struct part *p = &parts[-gp->id_or_neg_offset];
    radius = angular_smoothing_scale(x_cross, p->h);
  } else {
    radius = 0.0;
  }

  /* Loop over shells to update */
  for(int shell_nr=props->shell_nr_min; shell_nr<=props->shell_nr_max; shell_nr+=1) {
    if(a_cross > props->shell[shell_nr].amin && a_cross <= props->shell[shell_nr].amax) {
  
      /* Make sure this shell is available for updating */
      if(props->shell[shell_nr].state == shell_uninitialized)
        error("Attempt to update shell which has not been allocated");
      if(props->shell[shell_nr].state == shell_complete)
        error("Attempt to update shell which has been written out");

      /* Allocate storage for updates and set particle coordinates and radius */
      double *data = malloc(part_type_info->buffer_element_size);
      data[0] = theta;
      data[1] = phi;
      data[2] = radius;

      /* Loop over healpix maps which this particle type contributes to and find values to add */
      for(int i=0; i<part_type_info->nr_maps; i+=1) {
        int map_nr = part_type_info->map_index[i];
        data[3+i] = props->map_type[map_nr].update_map(e, props, gp, a_cross, x_cross);

#ifdef LIGHTCONE_MAP_CHECK_TOTAL
        /* Accumulate total quantity added to each map for consistency check */
        atomic_add_d(&props->shell[shell_nr].map[map_nr].total, data[3+i]);
#endif

      }

      /* Buffer the updates */
      particle_buffer_append(&(props->shell[shell_nr].buffer[gp->type]), data);

      /* Free update info */
      free(data);

    }
  } /* Next shell */
}


/**
 * @brief Compute memory used by lightcones on this rank
 *
 * @param props The #lightcone_props structure
 *
 */
void lightcone_report_memory_use(struct lightcone_props *props) {
  
  long long memuse_local[3];
  for(int i=0; i<3; i+=1)
    memuse_local[i] = 0;

  /* Accumulate memory used by particle buffers - one buffer per particle type */
  for(int i=0; i<swift_type_count; i+=1) {
    if(props->use_type[i])
      memuse_local[0] += particle_buffer_memory_use(props->buffer+i);
  }

  /* Accumulate memory used by map update buffers and pixel data */
  const int nr_maps = props->nr_maps;
  const int nr_shells = props->nr_shells;
  for(int shell_nr=0; shell_nr<nr_shells; shell_nr+=1) {

    /* Healpix map updates - one buffer per particle type per shell */
    for(int ptype=0; ptype<swift_type_count; ptype+=1) {
      memuse_local[1] += particle_buffer_memory_use(&(props->shell[shell_nr].buffer[ptype]));
    }

    /* Pixel data - one buffer per map per shell */
    for(int map_nr=0; map_nr<nr_maps; map_nr+=1) {
      struct lightcone_map *map = &(props->shell[shell_nr].map[map_nr]);
      if(map->data)memuse_local[2] += map->local_nr_pix*sizeof(double);
    }
  }

  /* Find min and max memory use over all nodes */
#ifdef WITH_MPI
  long long memuse_min[3];
  MPI_Reduce(memuse_local, memuse_min, 3, MPI_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
  long long memuse_max[3];
  MPI_Reduce(memuse_local, memuse_max, 3, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
  if(engine_rank==0) {
    message("lightcone %d: particle buffer bytes:   min=%lld, max=%lld", props->index, memuse_min[0], memuse_max[0]);
    message("lightcone %d: map update buffer bytes: min=%lld, max=%lld", props->index, memuse_min[1], memuse_max[1]);
    message("lightcone %d: map pixel data bytes:    min=%lld, max=%lld", props->index, memuse_min[2], memuse_max[2]);
  }
#else
    message("lightcone %d: particle buffer bytes:   %lld", props->index, memuse_local[0]);
    message("lightcone %d: map update buffer bytes: %lld", props->index, memuse_local[1]);
    message("lightcone %d: map pixel data bytes:    %lld", props->index, memuse_local[2]);
#endif    

}


/**
 * @brief Write out number of files per rank for this lightcone
 *
 * @param props The #lightcone_props structure
 *
 */
void lightcone_write_index(struct lightcone_props *props) {

  int comm_size = 1;
#ifdef WITH_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
#endif

  /* Collect current file index on each rank */
  int *current_file_on_rank = malloc(sizeof(int)*comm_size);
#ifdef WITH_MPI
  MPI_Gather(&props->current_file, 1, MPI_INT,
             current_file_on_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
#else
  current_file_on_rank[0] = props->current_file;
#endif

  if(engine_rank == 0) {

    /* Get the name of the index file */
    char fname[FILENAME_BUFFER_SIZE];
    int len = snprintf(fname, FILENAME_BUFFER_SIZE, "%s_index.hdf5", props->basename);
    if((len < 0) || (len >= FILENAME_BUFFER_SIZE))error("Failed to generate lightcone index filename");

    /* Create the file */
    hid_t file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    /* Write number of MPI ranks and number of files */
    hid_t group_id = H5Gcreate2(file_id, "Lightcone", H5P_DEFAULT,
                                H5P_DEFAULT, H5P_DEFAULT);
    io_write_attribute_i(group_id, "nr_mpi_ranks", comm_size);
    io_write_attribute(group_id, "final_file_on_rank", INT,
                       current_file_on_rank, comm_size);

    /* Write observer position and redshift limits */
    io_write_attribute(group_id, "observer_position", DOUBLE,
                       props->observer_position, 3);
    io_write_attribute_d(group_id, "minimum_redshift", props->z_min_for_particles);
    io_write_attribute_d(group_id, "maximum_redshift", props->z_max_for_particles);

    /* Write the number of shells and their radii */
    const int nr_shells = props->nr_shells;
    io_write_attribute_i(group_id, "nr_shells", nr_shells);
    double *shell_inner_radii = malloc(sizeof(double)*nr_shells);
    double *shell_outer_radii = malloc(sizeof(double)*nr_shells);
    for(int i=0; i<nr_shells; i+=1) {
      shell_inner_radii[i] = props->shell[i].rmin;
      shell_outer_radii[i] = props->shell[i].rmax;
    }
    io_write_attribute(group_id, "shell_inner_radii", DOUBLE,
                       shell_inner_radii, nr_shells);
    io_write_attribute(group_id, "shell_outer_radii", DOUBLE,
                       shell_outer_radii, nr_shells);
    free(shell_outer_radii);
    free(shell_inner_radii);

    H5Gclose(group_id);
    H5Fclose(file_id);
  }

  free(current_file_on_rank);
}
