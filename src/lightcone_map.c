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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* This object's header. */
#include "lightcone_map.h"

/* Local headers */
#include "align.h"
#include "common_io.h"
#include "error.h"
#include "exchange_structs.h"
#include "memuse.h"
#include "particle_buffer.h"
#include "restart.h"
#include "healpix_smoothing.h"

/* HDF5 */
#ifdef HAVE_HDF5
#include <hdf5.h>
#endif


#ifdef WITH_MPI
static int pixel_to_rank(struct lightcone_map *map, size_t pixel) {
  int rank = pixel / map->pix_per_rank;
  if(rank >= map->comm_size)rank = map->comm_size-1;
  return rank;
}
#endif

void lightcone_map_init(struct lightcone_map *map, const int nside,
                        const double r_min, const double r_max,
                        const size_t elements_per_block,
                        enum unit_conversion_factor units) {
  
  int comm_rank = 0, comm_size = 1;
#ifdef WITH_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif
  map->comm_size = comm_size;
  map->comm_rank = comm_rank;

  /* Initialise C++ code for smoothing on the sphere */
  map->smoothing_info = healpix_smoothing_init(nside);

  /* Initialise the data buffer for this map */
  particle_buffer_init(&map->buffer, sizeof(struct lightcone_map_contribution),
                       elements_per_block, "lightcone_map");

  /* Determine number of pixels in the map */  
  map->total_nr_pix = healpix_smoothing_get_npix(map->smoothing_info);

  /* Determine which pixels are stored on which rank:
     put pix_per_rank on each node with any extra assigned to
     the last node. This makes it easy to convert a pixel index
     to a node index. */
  map->pix_per_rank = map->total_nr_pix / comm_size;
  if(map->pix_per_rank==0)error("Must have healpix npix > number of MPI ranks!");
  if(comm_rank < comm_size-1)
    map->local_nr_pix = map->pix_per_rank;
  else
    map->local_nr_pix = map->total_nr_pix - (comm_size-1)*map->pix_per_rank;

  /* Store offset from local array index to global healpix pixel index  */
  map->local_pix_offset = map->pix_per_rank * comm_rank;
  
  /* Pixel data is initially not allocated */
  map->data = NULL;

  /* Record block size so we can re-initialise particle_buffer on restarting */
  map->elements_per_block = elements_per_block;
  
  /* Store resolution parameter, shell size, units */
  map->nside = nside;
  map->r_min = r_min;
  map->r_max = r_max;
  map->units = units;
}


/**
 * @brief Deallocate the lightcone_map pixel data
 *
 * @param map the #lightcone_map structure
 */
void lightcone_map_clean(struct lightcone_map *map) {
  
  healpix_smoothing_clean(map->smoothing_info);
  particle_buffer_free(&map->buffer);
  if(map->data)lightcone_map_free_pixels(map);
}


/**
 * @brief Allocate (and maybe initialize) the lightcone_map pixel data
 *
 * @param map the #lightcone_map structure
 * @param zero_pixels if true, set allocated pixels to zero
 */
void lightcone_map_allocate_pixels(struct lightcone_map *map, const int zero_pixels) {
  
  if(swift_memalign("lightcone_map_pixels", (void **) &map->data,
                    SWIFT_STRUCT_ALIGNMENT, sizeof(double)*map->local_nr_pix) != 0)
    error("Failed to allocate lightcone map pixel data");

  if(zero_pixels) {
    for(size_t i=0; i<map->local_nr_pix; i+=1)
      map->data[i] = 0.0;
  }

}


void lightcone_map_free_pixels(struct lightcone_map *map) {
  
  swift_free("lightcone_map_pixels", (void *) map->data);
  map->data = NULL;

}


/**
 * @brief Dump lightcone_map struct to the output stream.
 *
 * @param map the #lightcone_map structure
 * @param stream The stream to write to.
 */
void lightcone_map_struct_dump(const struct lightcone_map *map, FILE *stream) {

  /* Don't write the particle_buffer (must flush before dumping) */
  struct lightcone_map tmp = *map;
  memset(&tmp.buffer, 0, sizeof(struct particle_buffer));

  /* Don't write pointer to C++ smoothing info */
  tmp.smoothing_info = NULL;

  /* Write the struct */
  restart_write_blocks((void *) &tmp, sizeof(struct lightcone_map), 1, stream,
                       "lightcone_map", "lightcone_map");

  /* Write the pixel data if it is allocated */
  if(tmp.data)
    restart_write_blocks((void *) tmp.data, sizeof(double), tmp.local_nr_pix, 
                         stream, "lightcone_map_data", "lightcone_map_data");
}


/**
 * @brief Restore lightcone_map struct from the input stream.
 *
 * @param map the #lightcone_map structure
 * @param stream The stream to read from.
 */
void lightcone_map_struct_restore(struct lightcone_map *map, FILE *stream) {

  /* Read the struct */
  restart_read_blocks((void *)map, sizeof(struct lightcone_map), 1, stream,
                      NULL, "lightcone_map");

  /* Initialise the buffer for this map */
  particle_buffer_init(&map->buffer, sizeof(struct lightcone_map_contribution),
                       map->elements_per_block, "lightcone_map");
  
  /* Read the pixel data if it was allocated.
     map->data from the restart file is not a valid pointer now but we can
     check if it is not null to see if the pixel data block was written out. */
  if(map->data) {
    lightcone_map_allocate_pixels(map, /* zero_pixels = */ 0);
    restart_read_blocks((void *)map->data, sizeof(double), map->local_nr_pix,
                        stream, NULL, "lightcone_map");
  }

  /* Re-initialise C++ smoothing info */
  map->smoothing_info = healpix_smoothing_init(map->nside);

}


/**
 * @brief Mapper function for updating the healpix map
 *
 * @param map_data Pointer to an array of #lightcone_map_contribution
 * @param num_elements Number of elements in map_data
 * @param extra_data Pointer to the #lightcone_map struct
 *
 */
void healpix_smoothing_mapper(void *map_data, int num_elements,
                              void *extra_data) {
  
  /* Find the array of updates to apply */
  const struct lightcone_map_contribution *contr = 
    (struct lightcone_map_contribution *) map_data;

  /* Get a pointer to the lightcone_map struct to update */
  struct lightcone_map *map = (struct lightcone_map *) extra_data;

  /* Apply the updates */
  for(size_t i=0; i<num_elements; i+=1) {
    healpix_smoothing_add_to_map(map->smoothing_info, contr[i].pos,
                                 contr[i].radius, contr[i].value,
                                 map->local_pix_offset, map->local_nr_pix,
                                 map->data);
  }
}


struct buffer_count_data {
  size_t *count;
  size_t *offset;
  struct lightcone_map *map;
  struct lightcone_map_contribution *data;
};


#ifdef WITH_MPI
static void count_elements_to_send(void *map_data, int num_elements,
                                   void *extra_data) {
  
  /* Unpack extra data */
  struct buffer_count_data *send_data = (struct buffer_count_data *) extra_data;
  struct lightcone_map *map = send_data->map;

  /* Array of pointers to buffer blocks to process on this call */
  struct particle_buffer_block **blocks = (struct particle_buffer_block **) map_data;

  /* Loop over blocks to process */
  for(int block_nr=0; block_nr<num_elements; block_nr+=1) {

    /* Data in this block is an array of lightcone_map_contribution */
    struct lightcone_map_contribution *contr =
      (struct lightcone_map_contribution *) blocks[block_nr]->data;
    for(size_t i=0; i<blocks[block_nr]->num_elements; i+=1) {
  
      /* Determine which ranks this contribution goes to */
      size_t first_pixel, last_pixel;
      healpix_smoothing_get_pixel_range(map->smoothing_info, contr[i].pos,
                                        contr[i].radius, &first_pixel, &last_pixel);
      int first_dest = pixel_to_rank(map, first_pixel);
      int last_dest = pixel_to_rank(map, last_pixel);
      /* Increment the count for these destinations */
      for(int dest=first_dest; dest<=last_dest; dest+=1)      
        atomic_add(&(send_data->count[dest]), 1);
    }
  }
}


static void store_elements_to_send(void *map_data, int num_elements,
                                   void *extra_data) {
  
  /* Unpack extra data */
  struct buffer_count_data *send_data = (struct buffer_count_data *) extra_data;
  struct lightcone_map *map = send_data->map;

  /* Array of buffer blocks to process on this call */
  struct particle_buffer_block **blocks = (struct particle_buffer_block **) map_data;

  /* Loop over blocks to process */
  for(int block_nr=0; block_nr<num_elements; block_nr+=1) {

    /* Data in this block is an array of lightcone_map_contribution */
    struct lightcone_map_contribution *contr =
      (struct lightcone_map_contribution *) blocks[block_nr]->data;
    for(size_t i=0; i<blocks[block_nr]->num_elements; i+=1) {

      /* Determine which ranks this contribution goes to */
      size_t first_pixel, last_pixel;
      healpix_smoothing_get_pixel_range(map->smoothing_info, contr[i].pos,
                                        contr[i].radius, &first_pixel, &last_pixel);
      int first_dest = pixel_to_rank(map, first_pixel);
      int last_dest = pixel_to_rank(map, last_pixel);

      /* Loop over destination ranks */
      for(int dest=first_dest; dest<=last_dest; dest+=1) {
        /* Reserve a location to write this contribution to */
        size_t index = atomic_add(&(send_data->count[dest]), 1);
        /* Store the contribution to the send buffer */
        struct lightcone_map_contribution *from = contr+i;
        struct lightcone_map_contribution *to = send_data->data+send_data->offset[dest]+index;
        memcpy(to, from, sizeof(struct lightcone_map_contribution));
      }
    }
  }
}
#endif


/**
 * @brief Apply buffered updates to the healpix map
 *
 * @param map the #lightcone_map structure
 */
void lightcone_map_update_from_buffer(struct lightcone_map *map,
                                      struct threadpool *tp,
                                      const int verbose) {
  
#ifdef WITH_MPI

  /* Get MPI rank, number of ranks */
  const int comm_rank = map->comm_rank;
  const int comm_size = map->comm_size;

  /* Data to pass to threadpool map function */
  struct buffer_count_data send_data;
  send_data.map = map;

  /* Array to store number of updates to send to each rank*/
  send_data.count = malloc(comm_size*sizeof(size_t));
  for(int i=0; i<comm_size; i+=1)
    send_data.count[i] = 0;
  send_data.offset = NULL;
  send_data.data = NULL;

  /* Count elements to send to each rank */
  particle_buffer_threadpool_map(&map->buffer, tp, count_elements_to_send, &send_data);
    
  /* Find total number of updates */
  size_t total_nr_send = particle_buffer_num_elements(&map->buffer);
  if(verbose) {
    long long num_updates_local = total_nr_send;
    long long num_updates_total;
    MPI_Reduce(&num_updates_local, &num_updates_total, 1, MPI_LONG_LONG,
               MPI_SUM, 0, MPI_COMM_WORLD);
    if(comm_rank==0)
      message("total lightcone map updates to apply: %lld", num_updates_total);
  }

  /* Allocate send buffer */
  struct lightcone_map_contribution *sendbuf = 
    malloc(sizeof(struct lightcone_map_contribution)*total_nr_send);

  /* Compute offsets into send buffer */
  send_data.offset = malloc(comm_size*sizeof(size_t));
  send_data.offset[0] = 0;
  for(int i=1; i<comm_size; i+=1)
    send_data.offset[i] = send_data.offset[i-1] + send_data.count[i-1];

  /* Populate the send buffer */
  for(int i=0; i<comm_size; i+=1)
    send_data.count[i] = 0;
  send_data.data = sendbuf;
  particle_buffer_threadpool_map(&map->buffer, tp, store_elements_to_send, &send_data);

  /* Empty the particle buffer now that we copied the data from it */
  particle_buffer_empty(&map->buffer);

  /* Determine number of elements to receive */
  size_t *recv_count = malloc(comm_size*sizeof(size_t));
  MPI_Alltoall(send_data.count, sizeof(size_t), MPI_BYTE, recv_count, sizeof(size_t),
               MPI_BYTE, MPI_COMM_WORLD);
  size_t total_nr_recv = 0;
  for(int i=0; i<comm_size; i+=1)
    total_nr_recv += recv_count[i];
  
  /* Allocate receive buffer */
  struct lightcone_map_contribution *recvbuf = 
    malloc(sizeof(struct lightcone_map_contribution)*total_nr_recv);
  
  /* Exchange data */
  exchange_structs(send_data.count, sendbuf, recv_count, recvbuf,
                   sizeof(struct lightcone_map_contribution));

  /* Apply received updates to the healpix map */
  threadpool_map(tp, healpix_smoothing_mapper, recvbuf, total_nr_recv,
                 sizeof(struct lightcone_map_contribution), 1, map);

  /* Tidy up */
  free(send_data.count);
  free(send_data.offset);
  free(sendbuf);
  free(recv_count);
  free(recvbuf);

#else
  
  /* If not using MPI, we can update the map directly from the buffer */
  struct particle_buffer_block *block = NULL;
  struct lightcone_map_contribution *contr;
  size_t num_elements;
  do {
    particle_buffer_iterate(&map->buffer, &block, &num_elements, (void **) &contr);
    threadpool_map(tp, healpix_smoothing_mapper, contr, num_elements,
                   sizeof(struct lightcone_map_contribution), 
                   threadpool_auto_chunk_size, map);
  } while(block);
  particle_buffer_empty(&map->buffer);
    
#endif

}

#ifdef HAVE_HDF5
/**
 * @brief Write a lightcone map to a HDF5 file
 *
 * @param map the #lightcone_map structure
 * @param loc a HDF5 file or group identifier to write to
 * @param name the name of the dataset to create
 */
void lightcone_map_write(struct lightcone_map *map, const hid_t loc_id, const char *name,
                         const struct unit_system *internal_units,
                         const struct unit_system *snapshot_units) {

  /* Find unit conversion factor for this quantity */
  const double conversion_factor =
    units_conversion_factor(internal_units, snapshot_units, map->units);
  
  /* Convert units if necessary */
  if(conversion_factor != 1.0) {
    for(size_t i=0; i<map->local_nr_pix; i+=1)
      map->data[i] *= conversion_factor;
  }

  /* Create dataspace in memory corresponding to local pixels */
  const hsize_t mem_dims[1] = {(hsize_t) map->local_nr_pix};
  hid_t mem_space_id = H5Screate_simple(1, mem_dims, NULL);
  if(mem_space_id < 0)error("Unable to create memory dataspace");
  
  /* Create dataspace in the file corresponding to the full map */
  const hsize_t file_dims[1] = {(hsize_t) map->total_nr_pix};
  hid_t file_space_id = H5Screate_simple(1, file_dims, NULL);
  if(file_space_id < 0)error("Unable to create file dataspace");

  /* Select the part of the dataset in the file to write to */
#ifdef WITH_MPI
#ifdef HAVE_PARALLEL_HDF5
  const int comm_rank = map->comm_rank;
  const size_t pixel_offset = map->pix_per_rank * comm_rank;
  const hsize_t start[1] = {(hsize_t) pixel_offset};
  const hsize_t count[1] = {(hsize_t) map->local_nr_pix};
  if(H5Sselect_hyperslab(file_space_id, H5S_SELECT_SET, start, NULL, count, NULL) < 0)
    error("Unable to select part of file dataspace to write to");
#else
  error("Writing lightcone maps with MPI requires parallel HDF5");
#endif
#endif

  /* Create the dataset */
  hid_t dset_id = H5Dcreate(loc_id, name, H5T_NATIVE_DOUBLE, file_space_id,
                            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if(dset_id < 0)error("Unable to create dataset %s", name);
    
  /* Write attributes */
  io_write_attribute_i(dset_id, "nside", map->nside);
  io_write_attribute_l(dset_id, "number_of_pixels", map->total_nr_pix);
  io_write_attribute_s(dset_id, "pixel_ordering_scheme", "ring");
  io_write_attribute_d(dset_id, "comoving_inner_radius", map->r_min);
  io_write_attribute_d(dset_id, "comoving_outer_radius", map->r_max);

  /* Write unit conversion factors for this data set */
  char buffer[FIELD_BUFFER_SIZE] = {0};
  units_cgs_conversion_string(buffer, snapshot_units, map->units, 0.f);
  float baseUnitsExp[5];
  units_get_base_unit_exponents_array(baseUnitsExp, map->units);
  io_write_attribute_f(dset_id, "U_M exponent", baseUnitsExp[UNIT_MASS]);
  io_write_attribute_f(dset_id, "U_L exponent", baseUnitsExp[UNIT_LENGTH]);
  io_write_attribute_f(dset_id, "U_t exponent", baseUnitsExp[UNIT_TIME]);
  io_write_attribute_f(dset_id, "U_I exponent", baseUnitsExp[UNIT_CURRENT]);
  io_write_attribute_f(dset_id, "U_T exponent", baseUnitsExp[UNIT_TEMPERATURE]);
  io_write_attribute_f(dset_id, "h-scale exponent", 0.f);
  io_write_attribute_f(dset_id, "a-scale exponent", 0.f);
  io_write_attribute_s(dset_id, "Expression for physical CGS units", buffer);

  /* Write the actual number this conversion factor corresponds to */
  const double cgs_factor = units_cgs_conversion_factor(snapshot_units, map->units);
  io_write_attribute_d(dset_id,
                       "Conversion factor to CGS (not including cosmological corrections)",
                       cgs_factor);

  /* Set up property list for the write */
  hid_t h_plist_id = H5Pcreate(H5P_DATASET_XFER);
#if defined(WITH_MPI)
  if(H5Pset_dxpl_mpio(h_plist_id, H5FD_MPIO_COLLECTIVE) < 0)
    error("Unable to set collective transfer mode");
#endif

  /* Write the data */
  if(H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_space_id, file_space_id,
              h_plist_id, map->data) < 0)
    error("Unable to write dataset %s", name);

  /* Tidy up */
  H5Dclose(dset_id);
  H5Sclose(mem_space_id);
  H5Sclose(file_space_id);
  H5Pclose(h_plist_id);
  
}
#endif /* HAVE_HDF5*/
