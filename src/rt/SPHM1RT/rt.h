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
#ifndef SWIFT_RT_SPHM1RT_H
#define SWIFT_RT_SPHM1RT_H

#include "rt_properties.h"
#include "rt_struct.h"

#include <float.h>

/**
 * @file src/rt/SPHM1RT/rt.h
 * @brief Main header file for no radiative transfer scheme.
 * SPHM1RT method described in Chan+21: 2102.08404
 */

/**
 * @brief Returns the comoving radiation energy of a particle at the last
 * time the particle was kicked.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void
radiation_get_comoving_radiation_energy_multifrequency(
    const struct part* restrict p, float energy[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    energy[g] = p->rt_data.conserved[g].energy;
  }
}

/**
 * @brief Returns the physical radiation energy of a particle at the last
 * time the particle was kicked.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void
radiation_get_physical_radiation_energy_multifrequency(
    const struct part* restrict p, const struct cosmology* cosmo,
    float energy[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    energy[g] = p->rt_data.conserved[g].energy;
  }
}

/**
 * @brief Sets the comoving radiation energy of a particle
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
radiation_set_comoving_radiation_energy_multifrequency(
    struct part* p, const float energy[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.conserved[g].energy = energy[g];
  }
}

/**
 * @brief Sets the physical radiation energy of a particle
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
radiation_set_physical_radiation_energy_multifrequency(
    struct part* p, const struct cosmology* cosmo,
    const float energy[RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.conserved[g].energy = energy[g];
  }
}

/**
 * @brief Returns the comoving radiation flux of a particle at the last
 * time the particle was kicked.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void
radiation_get_comoving_radiation_flux_multifrequency(
    const struct part* restrict p, float fradtemp[3 * RT_NGROUPS]) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    fradtemp[0 + g * 3] = p->rt_data.conserved[g].flux[0];
    fradtemp[1 + g * 3] = p->rt_data.conserved[g].flux[1];
    fradtemp[2 + g * 3] = p->rt_data.conserved[g].flux[2];
  }
}

/**
 * @brief Returns the physical radiation flux of a particle at the last
 * time the particle was kicked.
 *
 * @param p Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static void
radiation_get_physical_radiation_flux_multifrequency(
    const struct part* restrict p, const struct cosmology* cosmo,
    float fradtemp[3 * RT_NGROUPS]) {

  for (int g = 0; g < RT_NGROUPS; g++) {
    fradtemp[0 + g * 3] = p->rt_data.conserved[g].flux[0];
    fradtemp[1 + g * 3] = p->rt_data.conserved[g].flux[1];
    fradtemp[2 + g * 3] = p->rt_data.conserved[g].flux[2];
  }
}

/**
 * @brief Sets the comoving radiation flux of a particle
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
radiation_set_comoving_radiation_flux_multifrequency(
    struct part* p, const float frad[3 * RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.conserved[g].flux[0] = frad[0 + g * 3];
    p->rt_data.conserved[g].flux[1] = frad[1 + g * 3];
    p->rt_data.conserved[g].flux[2] = frad[2 + g * 3];
  }
}

/**
 * @brief Sets the physical radiation flux of a particle
 *
 * @param p The particle of interest.
 * @param cosmo Cosmology data structure
 * @param u The physical internal energy
 */
__attribute__((always_inline)) INLINE static void
radiation_set_physical_radiation_flux_multifrequency(
    struct part* p, const struct cosmology* cosmo,
    const float frad[3 * RT_NGROUPS]) {
  for (int g = 0; g < RT_NGROUPS; g++) {
    p->rt_data.conserved[g].flux[0] = frad[0 + g * 3];
    p->rt_data.conserved[g].flux[1] = frad[1 + g * 3];
    p->rt_data.conserved[g].flux[2] = frad[2 + g * 3];
  }
}

/**
 * @brief Initialisation of the RT density loop related particle data.
 * Note: during initalisation (space_init), rt_reset_part and rt_init_part
 * are both called individually.
 */
__attribute__((always_inline)) INLINE static void rt_init_part(
    struct part* restrict p) {

  struct rt_part_data* rpd = &p->rt_data;

  float fradmag, fox;

  /* To avoid radiation reaching other dimension and violating conservation */
  for (int g = 0; g < RT_NGROUPS; g++) {
    if (hydro_dimension < 1.001f) {
      rpd->conserved[g].flux[1] = 0.0f;
    }
    if (hydro_dimension < 2.001f) {
      rpd->conserved[g].flux[2] = 0.0f;
    }
  }

  for (int g = 0; g < RT_NGROUPS; g++) {
    /* TK: avoid the radiation flux to violate causality. Impose a limit: F<Ec
     */
    fradmag =
        sqrtf(rpd->conserved[g].flux[0] * rpd->conserved[g].flux[0] +
              rpd->conserved[g].flux[1] * rpd->conserved[g].flux[1] +
              rpd->conserved[g].flux[2] * rpd->conserved[g].flux[2] + FLT_MIN);
    fox = max(1.f, fabsf(fradmag / (rpd->conserved[g].energy + FLT_MIN) /
                         (rpd->params.cred + FLT_MIN)));
    rpd->conserved[g].flux[0] = rpd->conserved[g].flux[0] / fox;
    rpd->conserved[g].flux[1] = rpd->conserved[g].flux[1] / fox;
    rpd->conserved[g].flux[2] = rpd->conserved[g].flux[2] / fox;
  }

  for (int g = 0; g < RT_NGROUPS; g++) {
    rpd->viscosity[g].divf = 0.0f;
    rpd->diffusion[g].graduradc[0] = 0.0f;
    rpd->diffusion[g].graduradc[1] = 0.0f;
    rpd->diffusion[g].graduradc[2] = 0.0f;
  }

  /* Some smoothing length multiples. */
  const float rho = hydro_get_comoving_density(p);
  const float rho_inv = 1.0f / rho; /* 1 / rho */

  /* Compute the "grad h" term */
  float rho_dh = p->density.rho_dh;

  const float omega_inv =
      1.f / (1.f + hydro_dimension_inv * p->h * rho_dh * rho_inv);

  /* Update variables. */
  rpd->force.f = omega_inv;
}

/**
 * @brief Reset of the RT hydro particle data not related to the density.
 * Note: during initalisation (space_init), rt_reset_part and rt_init_part
 * are both called individually. Also, if debugging checks are active, an
 * extra call to rt_reset_part is made in space_convert_rt_quantities() after
 * the zeroth time step is finished.
 */
__attribute__((always_inline)) INLINE static void rt_reset_part(
    struct part* restrict p) {

  struct rt_part_data* rpd = &p->rt_data;

  for (int g = 0; g < RT_NGROUPS; g++) {
    rpd->cdt[g].energy = 0.0f;
    rpd->cdt[g].flux[0] = 0.0f;
    rpd->cdt[g].flux[1] = 0.0f;
    rpd->cdt[g].flux[2] = 0.0f;
  }
}

/**
 * @brief First initialisation of the RT hydro particle data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_part(
    struct part* restrict p, const struct rt_props* restrict rt_props) {

  struct rt_part_data* rpd = &p->rt_data;

  for (int g = 0; g < RT_NGROUPS; g++) {
    rpd->viscosity[g].alpha = 1.0f;
    rpd->diffusion[g].alpha = 1.0f;
    rpd->params.chi[g] = rt_props->chi[g];
    rpd->viscosity[g].divf_previous_step = FLT_MIN;
  }

  /* We can get parameters for diffusion (force loop) */

  rpd->params.cred = rt_props->cred;

  rpd->force.f = 1.0f;

  rpd->dt = 1.0f;

  rt_init_part(p);
  rt_reset_part(p);
}

/**
 * @brief Initialisation of the RT density loop related star particle data.
 * Note: during initalisation (space_init), rt_reset_spart and rt_init_spart
 * are both called individually.
 */
__attribute__((always_inline)) INLINE static void rt_init_spart(
    struct spart* restrict sp) {}

/**
 * @brief Reset of the RT star particle data not related to the density.
 * Note: during initalisation (space_init), rt_reset_spart and rt_init_spart
 * are both called individually. Also, if debugging checks are active, an
 * extra call to rt_reset_spart is made in space_convert_rt_quantities() after
 * the zeroth time step is finished.
 */
__attribute__((always_inline)) INLINE static void rt_reset_spart(
    struct spart* restrict sp) {}

/**
 * @brief First initialisation of the RT star particle data.
 */
__attribute__((always_inline)) INLINE static void rt_first_init_spart(
    struct spart* restrict sp) {}

/**
 * @brief Split the RT data of a particle into n pieces
 *
 * @param p The #part.
 * @param n The number of pieces to split into.
 */
__attribute__((always_inline)) INLINE static void rt_split_part(struct part* p,
                                                                double n) {
  error("RT can't run with split particles for now.");
}

/**
 * @brief Exception handle a hydro part not having any neighbours in ghost task
 *
 * @param p The #part.
 */
__attribute__((always_inline)) INLINE static void rt_part_has_no_neighbours(
    struct part* p) {
  message("WARNING: found particle without neighbours");
};

/**
 * @brief Exception handle a star part not having any neighbours in ghost task
 *
 * @param p The #part.
 */
__attribute__((always_inline)) INLINE static void rt_spart_has_no_neighbours(
    struct spart* sp){};

/**
 * @brief Do checks/conversions on particles on startup.
 *
 * @param p The particle to work on
 * @param rtp The RT properties struct
 */
__attribute__((always_inline)) INLINE static void rt_convert_quantities(
    struct part* p, const struct rt_props* rtp){};

/**
 * @brief Computes the next radiative transfer time step size
 * of a given particle (during timestep tasks)
 *
 * @param p particle to work on
 * @param rt_props the RT properties struct
 * @param cosmo the cosmology
 */
__attribute__((always_inline)) INLINE static float rt_compute_timestep(
    const struct part* restrict p, const struct rt_props* restrict rt_props,
    const struct cosmology* restrict cosmo) {
  float dt = p->h * cosmo->a / (p->rt_data.params.cred + FLT_MIN) *
             rt_props->CFL_condition;

  return dt;
}

/**
 * @brief Compute the time-step length for an RT step of a particle.
 *
 * @param ti_beg Start of the time-step (on the integer time-line).
 * @param ti_end End of the time-step (on the integer time-line).
 * @param time_base Minimal time-step size on the time-line.
 * @param with_cosmology Are we running with cosmology integration?
 * @param cosmo The #cosmology object.
 *
 * @return The time-step size for the rt integration. (internal units).
 */
__attribute__((always_inline)) INLINE static double rt_part_dt(
    const integertime_t ti_beg, const integertime_t ti_end,
    const double time_base, const int with_cosmology,
    const struct cosmology* cosmo) {
  if (with_cosmology) {
    error("SPHM1RT with cosmology not implemented yet! :(");
    return 0.f;
  } else {
    return (ti_end - ti_beg) * time_base;
  }
}

/**
 * @brief Update the photon number of a particle, i.e. compute
 *  E^{n+1} = E^n + dt * dE_* / dt. This function finalises
 *  the injection step.
 *
 * @param p particle to work on
 * @param props struct #rt_props that contains global RT properties
 */
__attribute__((always_inline)) INLINE static void
rt_injection_update_photon_density(struct part* restrict p,
                                   struct rt_props* props) {}

/**
 * @brief Compute the photon emission rates for this stellar particle
 *        This function is called every time the spart is being reset
 *        (during start-up and during stars ghost if spart is active)
 *        and assumes that the photon emission rate is an intrinsic
 *        stellar property, i.e. doesn't depend on the environment.
 *
 * @param sp star particle to work on
 * @param time current system time
 * @param star_age age of the star *at the end of the step*
 * @param dt star time step
 * @param rt_props RT properties struct
 * @param phys_const physical constants struct
 * @param internal_units struct holding internal units
 */
__attribute__((always_inline)) INLINE static void
rt_compute_stellar_emission_rate(struct spart* restrict sp, double time,
                                 double star_age, double dt,
                                 const struct rt_props* rt_props,
                                 const struct phys_const* phys_const,
                                 const struct unit_system* internal_units) {}

/**
 * @brief finishes up the gradient computation
 *
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_end_gradient(
    struct part* restrict p) {
  struct rt_part_data* rpd = &p->rt_data;
  /* artificial diffusion for shock capturing */
  const float vsig_diss = rpd->params.cred;
  /* similar to Cullen & Dehnen 2010 switch */
  float divf, divf_previous_step, urad, viscosity_alpha, diffusion_alpha;
  float divf_dt, shockest, alphaflim, alpha_f_diss, alpha_f_diss_loc;
  float alpha_diss_loc, alpha_diss;

  for (int g = 0; g < RT_NGROUPS; g++) {
    divf = rpd->viscosity[g].divf;
    divf_previous_step = rpd->viscosity[g].divf_previous_step;
    urad = rpd->conserved[g].energy;
    viscosity_alpha = rpd->viscosity[g].alpha;
    diffusion_alpha = rpd->diffusion[g].alpha;
    divf_dt = (divf - divf_previous_step) / (rpd->dt + FLT_MIN);
    shockest = -p->h * p->h / (urad + FLT_MIN) / (vsig_diss + FLT_MIN) /
               (vsig_diss + FLT_MIN) * divf_dt * 200.f;
    alphaflim = max(shockest, 0.0f); /* should be positive or 0 */
    alpha_f_diss = viscosity_alpha;
    alpha_f_diss_loc = 0.0f;

    /* f diffusion only operates in compression */
    if (divf < 0.0f) {
      /* limit the diffusivity to Courant time step */
      alpha_f_diss_loc = min(alphaflim, 1.0f);
    }

    if (alpha_f_diss_loc > alpha_f_diss) {
      /* Reset the value of alpha to the appropriate value */
      alpha_f_diss = alpha_f_diss_loc;
    } else {
      /* Integrate the alpha forward in time to decay back to alpha = alpha_loc
       */
      alpha_f_diss =
          alpha_f_diss_loc +
          (alpha_f_diss - alpha_f_diss_loc) *
              expf(-rpd->dt * vsig_diss *
                   (1.f / (p->h + FLT_MIN) + rpd->params.chi[g] * p->rho));
    }

    /* alpha inspired by Price 2010: it should vanish where radiation energy
     * difference is small */
    alpha_diss_loc = 1.0f;
    alpha_diss = diffusion_alpha;
    if (alpha_diss_loc > alpha_diss) {
      /* Reset the value of alpha to the appropriate value */
      alpha_diss = alpha_diss_loc;
    } else {
      /* Integrate the alpha forward in time to decay back to alpha = alpha_loc
       */
      alpha_diss =
          alpha_diss_loc +
          (alpha_diss - alpha_diss_loc) *
              expf(-rpd->dt * vsig_diss *
                   (0.01f / (p->h + FLT_MIN) + rpd->params.chi[g] * p->rho));
    }

    /* Cap the dissipation to avoid instabilities */
    alpha_diss = min(alpha_diss, 1.0f);
    alpha_diss = max(alpha_diss, 0.0f);

    alpha_f_diss = min(alpha_f_diss, 1.0f);
    alpha_f_diss = max(alpha_f_diss, 0.0f);

    rpd->diffusion[g].alpha = alpha_diss;
    rpd->viscosity[g].alpha = alpha_f_diss;
  }
}

/**
 * @brief finishes up the transport step
 *
 * @param p particle to work on
 * @param dt the current time step of the particle
 */
__attribute__((always_inline)) INLINE static void rt_finalise_transport(
    struct part* restrict p, const double dt) {
  struct rt_part_data* rpd = &p->rt_data;

  for (int g = 0; g < RT_NGROUPS; g++) {
    rpd->conserved[g].energy += rpd->cdt[g].energy * dt;
    rpd->conserved[g].flux[0] += rpd->cdt[g].flux[0] * dt;
    rpd->conserved[g].flux[1] += rpd->cdt[g].flux[1] * dt;
    rpd->conserved[g].flux[2] += rpd->cdt[g].flux[2] * dt;
  }

  /* add frad source term implicitly */
  float dfrac;

  for (int g = 0; g < RT_NGROUPS; g++) {
    dfrac = -rpd->params.chi[g] * p->rho * rpd->params.cred;
    rpd->conserved[g].flux[0] *= expf(dfrac * dt);
    rpd->conserved[g].flux[1] *= expf(dfrac * dt);
    rpd->conserved[g].flux[2] *= expf(dfrac * dt);

    /* update urad */
    /* limiter to avoid negative urad */
    /* negative urad will make the dissipation (diffusion) unstable) */
    if (rpd->conserved[g].energy < 0.0f) {
      rpd->conserved[g].energy = FLT_MIN;
      rpd->conserved[g].flux[0] = FLT_MIN;
      rpd->conserved[g].flux[1] = FLT_MIN;
      rpd->conserved[g].flux[2] = FLT_MIN;
    }

    /* save next time step */
    rpd->viscosity[g].divf_previous_step = rpd->viscosity[g].divf;
  }

  rpd->dt = dt;

  /* To avoid radiation reaching other dimension and violating conservation */
  for (int g = 0; g < RT_NGROUPS; g++) {
    if (hydro_dimension < 1.001f) {
      rpd->conserved[g].flux[1] = 0.0f;
    }
    if (hydro_dimension < 2.001f) {
      rpd->conserved[g].flux[2] = 0.0f;
    }
  }
}

/**
 * @brief Do the thermochemistry on a particle.
 *
 * @param p particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_tchem(
    struct part* restrict p) {}

/**
 * @brief Clean the allocated memory inside the RT properties struct.
 *
 * @param props the #rt_props.
 */
__attribute__((always_inline)) INLINE static void rt_clean(
    struct rt_props* props) {}

#endif /* SWIFT_RT_SPHM1RT_H */
