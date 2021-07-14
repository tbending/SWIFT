#ifndef SWIFT_HEALPIX_SMOOTHING_H
#define SWIFT_HEALPIX_SMOOTHING_H

struct healpix_smoothing_info;

struct healpix_neighbour_info {
  size_t global_pix;
  double weight;
};
  
struct healpix_smoothing_info *healpix_smoothing_init(int nside, double gamma);

void healpix_smoothing_clean(struct healpix_smoothing_info *smooth_info);

size_t healpix_smoothing_get_npix(struct healpix_smoothing_info *smooth_info);

size_t healpix_smoothing_vec2pix(struct healpix_smoothing_info *smooth_info, 
                                 const double *pos);

void healpix_smoothing_vec2ang(struct healpix_smoothing_info *smooth_info, const double *pos,
                               double *theta, double *phi);

size_t healpix_smoothing_ang2pix(struct healpix_smoothing_info *smooth_info, 
                                 const double theta, const double phi);

double healpix_smoothing_get_max_pixrad(struct healpix_smoothing_info *smooth_info);

void healpix_smoothing_get_pixel_range(struct healpix_smoothing_info *smooth_info,
                                       const double theta, const double phi, const double radius,
                                       size_t *first_pixel, size_t *last_pixel);

void healpix_smoothing_find_neighbours(struct healpix_smoothing_info *smooth_info,
                                       const double theta, const double phi, const double radius,
                                       size_t *nr_ngb_out, struct healpix_neighbour_info **ngb);

#endif /* SWIFT_HEALPIX_SMOOTHING_H */
