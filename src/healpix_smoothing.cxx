#include <stdio.h>
#include <math.h>

#include "healpix_cxx/vec3.h"
#include "healpix_cxx/pointing.h"
#include "healpix_cxx/healpix_base.h"
#include "healpix_cxx/datatypes.h"

// This stores the healpix resolution info
static Healpix_Base2 hb;
static double max_pixrad;

// 2D Wendland C2 kernel (omitting normalisation)
// TODO: use same kernel as Swift?
static double projected_kernel(double r, double h) {
  
  const double w = pow(1.0-r/h, 4.0)*(1.0+4*r/h);
  if(w > 0)
    return w;
  else
    return 0.0;
}


extern "C" {

  void healpix_smoothing_init(int nside) {
    
    const Healpix_Ordering_Scheme scheme = RING;
    hb = Healpix_Base2(nside, scheme, SET_NSIDE);
    max_pixrad = hb.max_pixrad();

  }

  size_t healpix_smoothing_get_npix(void) {
    return (size_t) hb.Npix();
  }

  double healpix_smoothing_get_max_pixrad(void) {
    return hb.max_pixrad();
  }

  size_t healpix_smoothing_pixel_index(double *pos) {
    
    vec3 part_vec = vec3(pos[0], pos[1], pos[2]);
    return (size_t) hb.vec2pix(part_vec);
  }

  void healpix_smoothing_add_to_map(double *pos, double radius,
                                    double value, size_t local_pix_offset,
                                    size_t local_nr_pix, double *map_data) {
  
    // Get a normalized direction vector for this particle
    vec3 part_vec = vec3(pos[0], pos[1], pos[2]);
    part_vec.Normalize();

    // Small particles get added to a single pixel
    if(radius < max_pixrad) {
      int64 pixel = hb.vec2pix(part_vec);
      if((pixel >= local_pix_offset) && (pixel < local_pix_offset+local_nr_pix))
        map_data[pixel-local_pix_offset] += value;
      return;
    }

    // Find all pixels with centres within the angular radius
    // IMPORTANT: need to search a larger radius if kernel cutoff is > 1h
    std::vector<int64> pixels;
    hb.query_disc(pointing(part_vec), radius, pixels);
  
    // Loop over pixels within the radius
    double tot = 0.0;
    for(int64 pixel : pixels) {

      // Get direction vector to centre of this pixel
      vec3 pixel_vec = hb.pix2vec(pixel);
      pixel_vec.Normalize();

      // Find angle between this pixel centre and the particle
      const double angle = acos(dotprod(pixel_vec, part_vec));

      // Evaluate the kernel at this radius
      tot += projected_kernel(angle, radius);

    }
  
    // Now accumulate contributions to pixels
    for(int64 pixel : pixels) {
    
      // Get direction vector to centre of this pixel
      vec3 pixel_vec = hb.pix2vec(pixel);
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



