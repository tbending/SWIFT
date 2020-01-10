/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_LOGGER_MPI_HISTORY_H
#define SWIFT_LOGGER_MPI_HISTORY_H

#include "../config.h"

/* Standard includes */
#include <stdint.h>

/* Local include */
#include "part_type.h"
#include "error.h"


#if defined(WITH_LOGGER) && defined(WITH_MPI)

/* Forward declaration */
struct xpart;
struct part;
struct gpart;
struct spart;
struct bpart;

/**
 * @brief Contains the information concerning
 * a particle for the index files.
 */
struct logger_index_data {
  /* Id of the particle. */
  int64_t id;

  /* Offset of the particle in the file. */
  uint64_t offset;
};

/**
 * @brief Structure dealing with the particle
 * exchanged through MPI.
 */
struct logger_mpi_history {

  /* Number of elements currently stored */
  uint64_t size[swift_type_count];

  /* Size of the current buffer */
  size_t capacity[swift_type_count];

  /* Buffer containing the particles */
  struct logger_index_data *data[swift_type_count];
};

void logger_mpi_history_init(struct logger_mpi_history *hist,
                             int already_allocated);
void logger_mpi_history_clean(struct logger_mpi_history *hist);
void logger_mpi_history_log_part(struct logger_mpi_history *hist, struct part *p,
                                 struct xpart *xp);
void logger_mpi_history_log_spart(struct logger_mpi_history *hist, struct spart *sp);
void logger_mpi_history_log_gpart(struct logger_mpi_history *hist, struct gpart *gp);
void logger_mpi_history_log_bpart(struct logger_mpi_history *hist, struct bpart *bp);

#endif // WITH_LOGGER && WITH_MPI
#endif // SWIFT_LOGGER_MPI_HISTORY_H
