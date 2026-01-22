/**
 * GR-GSPH Pre-Interaction Implementation
 *
 * Extends SR-GSPH pre-interaction with metric-aware computations.
 */

#include "grgsph/gr_pre_interaction.hpp"
#include "srgsph/sr_primitive_recovery.hpp"
#include "simulation.hpp"
#include "parameters.hpp"
#include "periodic.hpp"
#include "openmp.hpp"
#include "kernel/kernel_function.hpp"
#include "bhtree.hpp"
#include "logger.hpp"


namespace sph {
namespace grgsph {

void GRPreInteraction::initialize(std::shared_ptr<SPHParameters> param)
{
    // Call base class initialization
    srgsph::PreInteraction::initialize(param);

    // Create metric if specified
    if (!m_metric && param->grgsph.metric_type != "none") {
        m_metric = create_metric(
            param->grgsph.metric_type,
            param->grgsph.bh_mass,
            param->grgsph.bh_spin
        );
    }
}

void GRPreInteraction::calculation(std::shared_ptr<Simulation> sim)
{
    // If no metric, fall back to SR behavior
    const MetricBase* metric = get_metric();
    if (!metric) {
        srgsph::PreInteraction::calculation(sim);
        return;
    }

    // Track how many times calculation has been called
    // The first 2 calls may happen before any force update, so we should
    // preserve initialized values
    static int call_count = 0;
    ++call_count;

    // GR-aware calculation
    if (m_first) {
        initial_smoothing(sim);
        m_first = false;
    }

    auto & particles = sim->get_particles();
    auto * periodic = sim->get_periodic().get();
    const int num = sim->get_particle_num();
    auto * kernel = sim->get_kernel().get();
    auto * tree = sim->get_tree().get();

    // MUSCL gradient arrays (if 2nd order enabled)
    std::vector<vec_t> * grad_p = nullptr;
    std::vector<vec_t> * grad_rho = nullptr;
    std::vector<vec_t> * grad_v[DIM] = {nullptr};

    if (m_is_2nd_order) {
        grad_p = &sim->get_vector_array("grad_pressure");
        grad_rho = &sim->get_vector_array("grad_density");
        grad_v[0] = &sim->get_vector_array("grad_velocity_0");
#if DIM >= 2
        grad_v[1] = &sim->get_vector_array("grad_velocity_1");
#endif
#if DIM == 3
        grad_v[2] = &sim->get_vector_array("grad_velocity_2");
#endif
    }

#pragma omp parallel for
    for (int i = 0; i < num; ++i) {
        auto & p_i = particles[i];

        std::vector<int> neighbor_list(m_neighbor_number * neighbor_list_size);

        const real search_r = p_i.sml * 6.0;
        const int n_neighbor_tmp = tree->neighbor_search(p_i, neighbor_list, particles, false);

        // Update smoothing length if iteration enabled
        if (m_iteration && !m_first_calculation) {
            p_i.sml = compute_smoothing_length(p_i, particles, neighbor_list,
                                               n_neighbor_tmp, periodic, kernel, p_i.sml);
        }

        // Ghost particles: only compute h
        if (p_i.is_ghost) {
            p_i.gradh = 1.0;
            continue;
        }

        // Compute particle volume and grad-h correction
        real gradh_i;
        const real Vp = compute_volume(p_i, particles, neighbor_list,
                                       n_neighbor_tmp, periodic, kernel, p_i.sml, gradh_i);
        p_i.gradh = gradh_i;

        // Number density
        const real N_new = p_i.nu / Vp;

        // === GR-AWARE LORENTZ FACTOR CALCULATION ===
        // Compute metric at particle position
        Metric31 metric_i;
        metric->compute(p_i.pos, metric_i);

        if (m_first_calculation) {
            // First timestep: compute conserved from initialized primitives
            // IMPORTANT: Use the initialized values (dens, pres, vel) as-is
            // Don't recompute them from volume - the initialization code set them correctly

            // Get coordinate velocity and transform to Eulerian for Lorentz factor
            real v_coord[3] = {0.0, 0.0, 0.0};
            for (int d = 0; d < DIM; ++d) {
                v_coord[d] = p_i.vel[d];
            }

            real V_eul[3];
            metric_i.coord_to_eulerian(v_coord, V_eul);

            // Compute proper Lorentz factor: Γ = 1/√(1 - γ_ij V^i V^j)
            const real gamma_lor = metric_i.lorentz_factor(V_eul);
            p_i.gamma_lor = gamma_lor;

            // Use INITIALIZED values for density and pressure (not recomputed)
            // p_i.dens was set by initialization as rest-frame density
            // p_i.pres was set by initialization
            const real rest_density = p_i.dens;
            const real u_internal = p_i.pres / ((m_gamma - 1.0) * std::max(rest_density, 1e-20));
            const real H = 1.0 + u_internal / (m_c_speed * m_c_speed) +
                           p_i.pres / (std::max(rest_density, 1e-20) * m_c_speed * m_c_speed);
            p_i.enthalpy = H;

            // Compute conserved variables with proper Lorentz factor
            // N = ρ_rest * Γ (lab-frame number density)
            p_i.N = rest_density * gamma_lor;
            p_i.S = p_i.vel * (gamma_lor * H);  // S^i = γH v^i

#if DIM == 1
            const real vel_t = p_i.vel_t;
            p_i.S_t = vel_t * (gamma_lor * H);
#else
            p_i.S_t = 0.0;
#endif

            const real X = m_gamma / (m_gamma - 1.0);
            p_i.e = (H * (X * gamma_lor * gamma_lor - 1.0) + 1.0) / (X * gamma_lor);

            // Keep initialized values
            p_i.ene = u_internal;

            // Sound speed (same formula as SR, uses proper thermodynamic h)
            p_i.sound = std::sqrt((m_gamma - 1.0) * (H - 1.0) / H) * m_c_speed;

        } else {
            // Normal timestep: recover primitives from conserved
            // BUT: on call_count <= 2, we're still in initialization phase
            // and should not do primitive recovery (no force update has happened yet)
            if (call_count <= 2) {
                // Skip primitive recovery on early calls - just update N
                // The conserved S and e were set correctly in first calculation
                // Keep primitives as-is, but update N from volume for next step
                p_i.N = N_new;
            } else {
                // Normal operation: recover primitives from conserved
                p_i.N = N_new;

                // Use standard primitive recovery (works in coordinate frame)
                // The recovered velocities are coordinate velocities v^i
                const auto prim = srgsph::PrimitiveRecovery::conserved_to_primitive_with_tangent(
                    p_i.S, p_i.S_t, p_i.e, p_i.N, m_gamma, m_c_speed);

                p_i.vel = prim.vel;
                p_i.vel_t = prim.vel_t;
                p_i.dens = prim.density;
                p_i.pres = prim.pressure;
                p_i.sound = prim.sound_speed;
                p_i.enthalpy = prim.enthalpy;

                // Recompute Lorentz factor with metric for force calculation
                real v_coord[3] = {0.0, 0.0, 0.0};
                for (int d = 0; d < DIM; ++d) {
                    v_coord[d] = prim.vel[d];
                }

                real V_eul[3];
                metric_i.coord_to_eulerian(v_coord, V_eul);
                p_i.gamma_lor = metric_i.lorentz_factor(V_eul);
            }
        }

        // Compute gradients for MUSCL reconstruction
        if (m_is_2nd_order) {
            vec_t grad_p_i(0.0);
            vec_t grad_rho_i(0.0);
            vec_t grad_vx_i(0.0), grad_vy_i(0.0), grad_vz_i(0.0);

            for (int n = 0; n < n_neighbor_tmp; ++n) {
                const int j = neighbor_list[n];
                if (i == j) continue;

                const auto & p_j = particles[j];
                const vec_t r_ij = periodic->calc_r_ij(p_i.pos, p_j.pos);
                const real r = std::abs(r_ij);

                if (r >= p_i.sml * 2.0) continue;

                const vec_t grad_W = kernel->dw(r_ij, r, p_i.sml);
                const real Vp_j = p_j.nu / p_j.N;

                grad_p_i += grad_W * (Vp_j * (p_j.pres - p_i.pres));
                grad_rho_i += grad_W * (Vp_j * (p_j.dens - p_i.dens));

                grad_vx_i += grad_W * (Vp_j * (p_j.vel[0] - p_i.vel[0]));
#if DIM >= 2
                grad_vy_i += grad_W * (Vp_j * (p_j.vel[1] - p_i.vel[1]));
#endif
#if DIM == 3
                grad_vz_i += grad_W * (Vp_j * (p_j.vel[2] - p_i.vel[2]));
#endif
            }

            (*grad_p)[i] = grad_p_i;
            (*grad_rho)[i] = grad_rho_i;
            (*grad_v[0])[i] = grad_vx_i;
#if DIM >= 2
            (*grad_v[1])[i] = grad_vy_i;
#endif
#if DIM == 3
            (*grad_v[2])[i] = grad_vz_i;
#endif
        }

        p_i.neighbor = n_neighbor_tmp;
    }

    m_first_calculation = false;
}

} // namespace grgsph
} // namespace sph
