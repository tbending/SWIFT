//extern "C" {

  void healpix_smoothing_init(int nside);

  size_t healpix_smoothing_get_npix(void);

  size_t healpix_smoothing_pixel_index(double *pos);

  double healpix_smoothing_get_max_pixrad(void);

  void healpix_smoothing_add_to_map(double *pos, double radius,
                                    double value, size_t local_pix_offset,
                                    size_t local_nr_pix, double *map_data);
//}
