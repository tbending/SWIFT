#include <stdio.h>
#include <math.h>

#include "healpix_cxx/vec3.h"
#include "healpix_cxx/pointing.h"
#include "healpix_cxx/healpix_base.h"
#include "healpix_cxx/datatypes.h"


// 2D Wendland C2 kernel (omitting normalisation)
// TODO: use same kernel as Swift?
static double projected_kernel(double r, double h) {
  
  const double w = pow(1.0-r/h, 4.0)*(1.0+4*r/h);
  if(w > 0)
    return w;
  else
    return 0.0;
}


struct healpix_smoothing_info {
  int nside;
  double max_pixrad;
  Healpix_Base2 healpix_base;
};


extern "C" {

  struct healpix_smoothing_info *healpix_smoothing_init(int nside) {
    
    const Healpix_Ordering_Scheme scheme = RING;
    healpix_smoothing_info *smooth_info = new struct healpix_smoothing_info;
    smooth_info->nside = nside;
    smooth_info->healpix_base = Healpix_Base2(nside, scheme, SET_NSIDE);
    smooth_info->max_pixrad = smooth_info->healpix_base.max_pixrad();
    return smooth_info;
  }

  void healpix_smoothing_clean(struct healpix_smoothing_info *smooth_info) {
    delete smooth_info;
  }

  size_t healpix_smoothing_get_npix(struct healpix_smoothing_info *smooth_info) {
    return (size_t) smooth_info->healpix_base.Npix();
  }

  double healpix_smoothing_get_max_pixrad(struct healpix_smoothing_info *smooth_info) {
    return smooth_info->healpix_base.max_pixrad();
  }

  size_t healpix_smoothing_pixel_index(struct healpix_smoothing_info *smooth_info, double *pos) {
    
    vec3 part_vec = vec3(pos[0], pos[1], pos[2]);
    return (size_t) smooth_info->healpix_base.vec2pix(part_vec);
  }

  void healpix_smoothing_add_to_map(struct healpix_smoothing_info *smooth_info,
                                    double *pos, double radius,
                                    double value, size_t local_pix_offset,
                                    size_t local_nr_pix, double *map_data) {
  
    // Get a normalized direction vector for this particle
    vec3 part_vec = vec3(pos[0], pos[1], pos[2]);
    part_vec.Normalize();

    // Small particles get added to a single pixel
    if(radius < smooth_info->max_pixrad) {
      int64 pixel = smooth_info->healpix_base.vec2pix(part_vec);
      if((pixel >= local_pix_offset) && (pixel < local_pix_offset+local_nr_pix))
        map_data[pixel-local_pix_offset] += value;
      return;
    }

    // Find all pixels with centres within the angular radius
    // IMPORTANT: need to search a larger radius if kernel cutoff is > 1h
    std::vector<int64> pixels;
    smooth_info->healpix_base.query_disc(pointing(part_vec), radius, pixels);
  
    // Loop over pixels within the radius
    double tot = 0.0;
    for(int64 pixel : pixels) {

      // Get direction vector to centre of this pixel
      vec3 pixel_vec = smooth_info->healpix_base.pix2vec(pixel);
      pixel_vec.Normalize();

      // Find angle between this pixel centre and the particle
      const double angle = acos(dotprod(pixel_vec, part_vec));

      // Evaluate the kernel at this radius
      tot += projected_kernel(angle, radius);

    }
  
    // Now accumulate contributions to pixels
    for(int64 pixel : pixels) {
    
      // Get direction vector to centre of this pixel
      vec3 pixel_vec = smooth_info->healpix_base.pix2vec(pixel);
      pixel_vec.Normalize();

      // Find angle between this pixel centre and the particle
      const double angle = acos(dotprod(pixel_vec, part_vec));

      // Evaluate the kernel at this radius
      const double weight = projected_kernel(angle, radius) / tot;

      // Add contribution to the local part of the map
      if((pixel >= local_pix_offset) && (pixel < local_pix_offset+local_nr_pix))
        map_data[pixel-local_pix_offset] += weight*value;
      
    }
  }

} // extern "C"



