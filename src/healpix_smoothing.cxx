#include <stdio.h>
#include <math.h>

#include "config.h"

#ifdef HAVE_HEALPIX_CXX
#include "healpix_cxx/vec3.h"
#include "healpix_cxx/pointing.h"
#include "healpix_cxx/healpix_base.h"
#include "healpix_cxx/datatypes.h"
#endif

extern "C" {
#include "atomic.h"
#include "healpix_smoothing.h"
#include "projected_kernel.h"
}


struct healpix_smoothing_info {
  int nside;
  double max_pixrad;
  double kernel_gamma;
#ifdef HAVE_HEALPIX_CXX
  Healpix_Base2 healpix_base;
#endif
  struct projected_kernel_table kernel;
};


extern "C" {

  struct healpix_smoothing_info *healpix_smoothing_init(int nside, double gamma) {
#ifdef HAVE_HEALPIX_CXX
    const Healpix_Ordering_Scheme scheme = RING;
    healpix_smoothing_info *smooth_info = new struct healpix_smoothing_info;
    smooth_info->nside = nside;
    smooth_info->healpix_base = Healpix_Base2(nside, scheme, SET_NSIDE);
    smooth_info->max_pixrad = smooth_info->healpix_base.max_pixrad();
    smooth_info->kernel_gamma = gamma;
    projected_kernel_init(&smooth_info->kernel);
    return smooth_info;
#else
    error("Code was compiled without Healpix C++ library");
    return NULL;
#endif
  }

  void healpix_smoothing_clean(struct healpix_smoothing_info *smooth_info) {
#ifdef HAVE_HEALPIX_CXX
    projected_kernel_clean(&smooth_info->kernel);
    delete smooth_info;
#else
    error("Code was compiled without Healpix C++ library");
#endif
  }

  size_t healpix_smoothing_get_npix(struct healpix_smoothing_info *smooth_info) {
#ifdef HAVE_HEALPIX_CXX
    return (size_t) smooth_info->healpix_base.Npix();
#else
    error("Code was compiled without Healpix C++ library");
    return 0;
#endif
  }

  double healpix_smoothing_get_max_pixrad(struct healpix_smoothing_info *smooth_info) {
#ifdef HAVE_HEALPIX_CXX
    return smooth_info->healpix_base.max_pixrad();
#else
    error("Code was compiled without Healpix C++ library");
    return -1.0;
#endif
  }

  size_t healpix_smoothing_vec2pix(struct healpix_smoothing_info *smooth_info, 
                                   const double *pos) {
#ifdef HAVE_HEALPIX_CXX
    vec3 part_vec = vec3(pos[0], pos[1], pos[2]);
    return (size_t) smooth_info->healpix_base.vec2pix(part_vec);
#else
    error("Code was compiled without Healpix C++ library");
    return 0;
#endif
  }

  void healpix_smoothing_vec2ang(struct healpix_smoothing_info *smooth_info,
                                 const double *pos, double *theta, double *phi) {
#ifdef HAVE_HEALPIX_CXX
    pointing p = pointing(vec3(pos[0], pos[1], pos[2]));
    *theta = p.theta;
    *phi = p.phi;
#else
    error("Code was compiled without Healpix C++ library");
#endif
  }

  size_t healpix_smoothing_ang2pix(struct healpix_smoothing_info *smooth_info,
                                   const double theta, const double phi) {
#ifdef HAVE_HEALPIX_CXX
    pointing p = pointing(theta, phi);
    return (size_t) smooth_info->healpix_base.ang2pix(p);
#else
    error("Code was compiled without Healpix C++ library");
    return 0;
#endif
  }

  void healpix_smoothing_get_pixel_range(struct healpix_smoothing_info *smooth_info,
                                         const double theta, const double phi, 
                                         const double radius, size_t *first_pixel,
                                         size_t *last_pixel) {
#ifdef HAVE_HEALPIX_CXX
    pointing p = pointing(theta, phi);

    // Input radius is the angle corresponding to the smoothing length,
    // so we need to find all pixels within kernel_gamma*radius
    const double search_radius = radius*smooth_info->kernel_gamma;

    // Small particles get added to a single pixel
    if(search_radius < smooth_info->max_pixrad) {
      int64 pixel = smooth_info->healpix_base.ang2pix(p);
      *first_pixel = pixel;
      *last_pixel = pixel;
      return;
    }

    // Find all pixels with centres within the angular radius
    std::vector<int64> pixels;
    smooth_info->healpix_base.query_disc(p, search_radius, pixels);
    *first_pixel = pixels[0];
    *last_pixel = pixels[0];
    for(size_t i=0; i < pixels.size(); i++) {
      if(pixels[i] < *first_pixel)*first_pixel = pixels[i];
      if(pixels[i] > *last_pixel)*last_pixel = pixels[i];
    }
#else
    error("Code was compiled without Healpix C++ library");
#endif
    return;
  }


  void healpix_smoothing_find_neighbours(struct healpix_smoothing_info *smooth_info,
                                         const double theta, const double phi, const double radius,
                                         size_t *nr_ngb_out, struct healpix_neighbour_info **ngb_out) {
#ifdef HAVE_HEALPIX_CXX
    pointing p = pointing(theta, phi);

    // Input radius is the angle corresponding to the smoothing length,
    // so we need to find all pixels within kernel_gamma*radius
    const double search_radius = radius*smooth_info->kernel_gamma;
    
    // Find all pixels with centres within the angular radius
    std::vector<int64> pixels;
    smooth_info->healpix_base.query_disc(p, search_radius, pixels);
    size_t nr_ngb = pixels.size();

    // Allocate output array
    struct healpix_neighbour_info *ngb = 
      (struct healpix_neighbour_info *) malloc(sizeof(struct healpix_neighbour_info)*nr_ngb);

    // Particle direction vector
    vec3 part_vec = p.to_vec3();

    // Compute weighting factors
    double tot = 0.0;
    for(size_t i=0; i<nr_ngb; i+=1) {

      // Get direction vector to centre of this pixel
      vec3 pixel_vec = smooth_info->healpix_base.pix2vec(pixels[i]);

      // Find angle between this pixel centre and the particle.
      // Dot product may be a tiny bit greater than one due to rounding error
      const double dp = dotprod(pixel_vec, part_vec);
      const double angle = dp < 1.0 ? acos(dp) : 0.0;
      
      // Evaluate the kernel at this radius
      ngb[i].global_pix = (size_t) pixels[i];
      ngb[i].weight = projected_kernel_eval(&smooth_info->kernel, angle/radius);
      tot += ngb[i].weight;
    }

    // Normalize the weights
    for(size_t i=0; i<nr_ngb; i+=1)
      ngb[i].weight /= tot;

    *nr_ngb_out = nr_ngb;
    *ngb_out = ngb;

#else
    error("Code was compiled without Healpix C++ library");
#endif
  }

} // extern "C"



