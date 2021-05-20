
struct healpix_smoothing_info;

struct healpix_smoothing_info *healpix_smoothing_init(int nside);

void healpix_smoothing_clean(struct healpix_smoothing_info *smooth_info);

size_t healpix_smoothing_get_npix(struct healpix_smoothing_info *smooth_info);

size_t healpix_smoothing_pixel_index(struct healpix_smoothing_info *smooth_info, double *pos);

double healpix_smoothing_get_max_pixrad(struct healpix_smoothing_info *smooth_info);

void healpix_smoothing_add_to_map(struct healpix_smoothing_info *smooth_info,
                                  const double *pos, const double radius,
                                  const double value, const size_t local_pix_offset,
                                  const size_t local_nr_pix, double *map_data);

