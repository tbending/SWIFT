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

#include <stdio.h>
#include <math.h>
#include <limits.h>

/* SWIFT configuration */
#include "config.h"

/* HEALPix C API */
#ifdef HAVE_CHEALPIX
#include "chealpix.h"
#endif

/* Local includes */
#include "atomic.h"
#include "chealpix_smoothing.h"
#include "projected_kernel.h"
#include "healpix_util.h"

void chealpix_smoothing_init(struct chealpix_smoothing_info *smooth_info,
                             int nside, double gamma) {

#ifdef HAVE_CHEALPIX
  /* Check on data types: the HEALPix C API uses longs to store pixel indexes, so check
     that this is sufficient for the specified nside. */
  long long npix = (12ll*nside)*nside;
  if(npix > LONG_MAX)error("Number of HEALPix pixels does not fit in a long");
  smooth_info->nside = nside;
  smooth_info->max_pixrad = max_pixrad(nside);
  smooth_info->kernel_gamma = gamma;
  projected_kernel_init(&smooth_info->kernel);
#else
  error("Code was compiled without Healpix C library");
  return NULL;
#endif
}

void chealpix_smoothing_clean(struct chealpix_smoothing_info *smooth_info) {
  projected_kernel_clean(&smooth_info->kernel);
}

size_t chealpix_smoothing_get_npix(struct chealpix_smoothing_info *smooth_info) {
  const int nside = smooth_info->nside;
  return (12ll*nside)*nside;
}

double chealpix_smoothing_get_max_pixrad(struct chealpix_smoothing_info *smooth_info) {
  return smooth_info->max_pixrad;
}

size_t chealpix_smoothing_vec2pix(struct chealpix_smoothing_info *smooth_info, 
                                  const double *pos) {
#ifdef HAVE_CHEALPIX
  long ipring;
  vec2pix_ring(smooth_info->nside, pos, &ipring);
  return (size_t) ipring;
#else
  error("Code was compiled without Healpix C library");
  return 0;
#endif
}

void chealpix_smoothing_vec2ang(struct chealpix_smoothing_info *smooth_info,
                                const double *pos, double *theta, double *phi) {
#ifdef HAVE_CHEALPIX
  vec2ang(pos, theta, phi);
#else
  error("Code was compiled without Healpix C library");
#endif
}

size_t chealpix_smoothing_ang2pix(struct chealpix_smoothing_info *smooth_info,
                                  const double theta, const double phi) {
#ifdef HAVE_CHEALPIX
  long ipring;
  ang2pix_ring(smooth_info->nside, theta, phi, &ipring);
  return (size_t) ipring;
#else
  error("Code was compiled without Healpix C library");
  return 0;
#endif
}

void chealpix_smoothing_get_pixel_range(struct chealpix_smoothing_info *smooth_info,
                                        const double theta, const double phi, 
                                        const double radius, size_t *first_pixel,
                                        size_t *last_pixel) {
#ifdef CHAVE_HEALPIX

  /* Convert spherical coords to a vector */
  double vec[3];
  ang2vec(theta, phi, vec);

  /* Input radius is the angle corresponding to the smoothing length,
     so we need to find all pixels within kernel_gamma*radius */
  const double search_radius = radius*smooth_info->kernel_gamma;

  /* Small particles get added to a single pixel */
  if(search_radius < smooth_info->max_pixrad) {
    *first_pixel = chealpix_smoothing_ang2pix(smooth_info, theta, phi);
    *last_pixel = *first_pixel;
    return;
  }

  /* Find all pixels with centres within the angular radius */
  long long pix_min;
  long long pix_max;
  query_disc_range(smooth_info->nside, vec, search_radius,
                   &pix_min, &pix_max, NULL, NULL);
  *first_pixel = pix_min;
  *last_pixel  = pix_max;
  return;

#else
  error("Code was compiled without Healpix C library");
  return;
#endif
}


void chealpix_smoothing_find_neighbours(struct chealpix_smoothing_info *smooth_info,
                                        const double theta, const double phi, const double radius,
                                        size_t *nr_ngb_out, struct chealpix_neighbour_info **ngb_out) {
#ifdef HAVE_CHEALPIX

  /* Convert spherical coords to a vector */
  double part_vec[3];
  ang2vec(theta, phi, part_vec);

  /* Input radius is the angle corresponding to the smoothing length,
     so we need to find all pixels within kernel_gamma*radius */
  const double search_radius = radius*smooth_info->kernel_gamma;
    
  /* Find all pixels with centres within the angular radius */
  int nr_ranges;
  struct pixel_range *range;
  long long pix_min;
  long long pix_max;
  query_disc_range(smooth_info->nside, part_vec, search_radius,
    &pix_min, &pix_max, &nr_ranges, &range);
  
  /* Find total number of neighbours */
  size_t nr_ngb = 0;
  for(int range_nr=0; range_nr<nr_ranges; range_nr+=1)
    nr_ngb += range[range_nr].last - range[range_nr].first + 1;

  /* Allocate output array */
  struct chealpix_neighbour_info *ngb =
    (struct chealpix_neighbour_info *) malloc(sizeof(struct chealpix_neighbour_info)*nr_ngb);

  /* Compute weighting factors */
  double tot = 0.0;
  size_t ngb_nr = 0;
  for(int range_nr=0; range_nr<nr_ranges; range_nr+=1) {
    for(size_t pix=range[range_nr].first; pix<=range[range_nr].last; pix+=1) {
      
      /* Get direction vector to centre of this pixel */
      double pixel_vec[3];
      pix2vec_ring(smooth_info->nside, (long) pix, pixel_vec);

      /* Find angle between this pixel centre and the particle.
         Dot product may be a tiny bit greater than one due to rounding error */
      const double dp = (pixel_vec[0]*part_vec[0] +
                         pixel_vec[1]*part_vec[1] +
                         pixel_vec[2]*part_vec[2]);
      const double angle = dp < 1.0 ? acos(dp) : 0.0;

      /* Evaluate the kernel at this radius */
      ngb[ngb_nr].global_pix = (size_t) pix;
      ngb[ngb_nr].weight = projected_kernel_eval(&smooth_info->kernel, angle/radius);
      tot += ngb[ngb_nr].weight;
      ngb_nr += 1;
      /* Next pixel in this range */
    }
    /* Next range of pixels */
  }

  /* Free the array of ranges from query_disc_range() */
  free(range);

  /* Normalize the weights */
  for(size_t i=0; i<nr_ngb; i+=1)
    ngb[i].weight /= tot;

  *nr_ngb_out = nr_ngb;
  *ngb_out = ngb;

#else
  error("Code was compiled without Healpix C library");
#endif
}




