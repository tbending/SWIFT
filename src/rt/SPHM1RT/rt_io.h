/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Tsang Keung Chan (chantsangkeung@gmail.com)
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_IO_SPHM1RT_H
#define SWIFT_RT_IO_SPHM1RT_H

/*#include "io_properties.h" */

/**
 * @file src/rt/SPHM1RT/rt_io.h
 * @brief Main header file for no radiative transfer scheme IO routines.
 * SPHM1RT method described in Chan+21: 2102.08404
 */

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int rt_read_particles(const struct part* parts,
                                    struct io_props* list) {

  /* List what we want to read */

  char fieldname[30];
  int count = 0;
  for (int phg = 0; phg < RT_NGROUPS; phg++) {
    sprintf(fieldname, "PhotonEnergiesGroup%d", phg + 1);
    list[count++] =
        io_make_input_field(fieldname, FLOAT, 1, OPTIONAL, UNIT_CONV_ENERGY,
                            parts, rt_data.conserved[phg].energy);
    sprintf(fieldname, "PhotonFluxesGroup%d", phg + 1);
    list[count++] = io_make_input_field(fieldname, FLOAT, 3, OPTIONAL,
                                        UNIT_CONV_RADIATION_FLUX, parts,
                                        rt_data.conserved[phg].flux);
  }

  return count;
}

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int rt_read_stars(const struct spart* sparts,
                                struct io_props* list) {
  return 0;
}

/**
 * @brief Extract photon energies of conserved struct for all photon groups
 */
INLINE static void rt_convert_conserved_photon_energies(
    const struct engine* engine, const struct part* part,
    const struct xpart* xpart, float* ret) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    ret[g] = part->rt_data.conserved[g].energy;
  }
}

/**
 * @brief Extract photon energies of conserved struct for all photon groups
 */
INLINE static void rt_convert_conserved_photon_fluxes(
    const struct engine* engine, const struct part* part,
    const struct xpart* xpart, float* ret) {

  int i = 0;
  for (int g = 0; g < RT_NGROUPS; g++) {
    ret[i++] = part->rt_data.conserved[g].flux[0];
    ret[i++] = part->rt_data.conserved[g].flux[1];
    ret[i++] = part->rt_data.conserved[g].flux[2];
  }
}

/**
 * @brief Creates additional output fields for the radiative
 * transfer data of hydro particles.
 */
INLINE static int rt_write_particles(const struct part* parts,
                                     struct io_props* list) {

  int num_elements = 2;

  list[0] = io_make_output_field_convert_part(
      "PhotonEnergies", FLOAT, RT_NGROUPS, UNIT_CONV_ENERGY, 0, parts,
      /*xparts=*/NULL, rt_convert_conserved_photon_energies,
      "Photon Energies (all groups)");
  list[1] = io_make_output_field_convert_part(
      "PhotonFluxes", FLOAT, 3 * RT_NGROUPS, UNIT_CONV_RADIATION_FLUX, 0, parts,
      /*xparts=*/NULL, rt_convert_conserved_photon_fluxes,
      "Photon Fluxes (all groups; x, y, and z coordinates)");

  return num_elements;
}

/**
 * @brief Creates additional output fields for the radiative
 * transfer data of star particles.
 */
INLINE static int rt_write_stars(const struct spart* sparts,
                                 struct io_props* list) {
  return 0;
}

/**
 * @brief Write the RT model properties to the snapshot.
 *
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 * @param e The engine
 * @param internal_units The internal unit system
 * @param snapshot_units Units used for the snapshot
 * @param rtp The #rt_props
 */
INLINE static void rt_write_flavour(hid_t h_grp, hid_t h_grp_columns,
                                    const struct engine* e,
                                    const struct unit_system* internal_units,
                                    const struct unit_system* snapshot_units,
                                    const struct rt_props* rtp) {

#if defined(HAVE_HDF5)

  /* Write scheme name */
  io_write_attribute_s(h_grp, "RT Scheme", RT_IMPLEMENTATION);
  
  /* Write photon group counts */
  /* ------------------------- */
  io_write_attribute_i(h_grp, "PhotonGroupNumber", RT_NGROUPS);

#endif /* HAVE_HDF5 */

}

#endif /* SWIFT_RT_IO_SPHM1RT_H */
