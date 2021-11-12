#ifndef SWIFT_CHEALPIX_SMOOTHING_H
#define SWIFT_CHEALPIX_SMOOTHING_H

/* Local includes */
#include "projected_kernel.h"

struct chealpix_smoothing_info {
  int nside;
  double max_pixrad;
  double kernel_gamma;
  struct projected_kernel_table kernel;
};

struct chealpix_neighbour_info {
  size_t global_pix;
  double weight;
};

void chealpix_smoothing_init(struct chealpix_smoothing_info *smooth_info,
                             int nside, double gamma);

void chealpix_smoothing_clean(struct chealpix_smoothing_info *smooth_info);

size_t chealpix_smoothing_get_npix(struct chealpix_smoothing_info *smooth_info);

size_t chealpix_smoothing_vec2pix(struct chealpix_smoothing_info *smooth_info,
                                  const double *pos);

void chealpix_smoothing_vec2ang(struct chealpix_smoothing_info *smooth_info,
                                const double *pos, double *theta, double *phi);

size_t chealpix_smoothing_ang2pix(struct chealpix_smoothing_info *smooth_info,
                                  const double theta, const double phi);

double chealpix_smoothing_get_max_pixrad(
    struct chealpix_smoothing_info *smooth_info);

void chealpix_smoothing_get_pixel_range(
    struct chealpix_smoothing_info *smooth_info, const double theta,
    const double phi, const double radius, size_t *first_pixel,
    size_t *last_pixel);

void chealpix_smoothing_find_neighbours(
    struct chealpix_smoothing_info *smooth_info, const double theta,
    const double phi, const double radius, size_t *nr_ngb_out,
    struct chealpix_neighbour_info **ngb);

#endif /* SWIFT_CHEALPIX_SMOOTHING_H */
