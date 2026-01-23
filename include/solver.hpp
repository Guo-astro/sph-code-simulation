#pragma once

#include <memory>
#include <unordered_map>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/any.hpp>

#include "defines.hpp"
#include "output_manager.hpp"
#include "output_metadata.hpp"
#include "units.hpp"

namespace sph
{

struct SPHParameters;
class Simulation;
class LaneEmdenRelaxation;
class PolytropicSlabRelaxation;
class PolytropicSlab2DRelaxation;
class KoyamaInutsukaRelaxation;
class IsothermalRelaxation;

class Module;

// Forward declaration for external forces
namespace external_forces {
    class PointMassBlackHole;
}

enum struct Sample {
    ShockTube,
    ShockTube2D,
    ShockTube3D,
    Vacuum,
    StrongShock,
    PressureEquilibrium,
    GreshoChanVortex,
    PairingInstability,
    HydroStatic,
    KHI,
    ISMCooling1D,
    Evrard,
    EvrardColdCollapse,
    LaneEmden,
    LaneEmdenCylinder,   // 3D cylindrical Lane-Emden (radial gravity in xy-plane)
    PolytropicSlab2D,    // 2D planar Lane-Emden slab (gravity in y-direction)
    Sedov,
    SRSod,
    SRTangentVelocity,
    SRRosswog,           // Rosswog (2010) SR-SPH benchmark tests
    GRSchwarzschildShock,// GR-GSPH Schwarzschild radial shock tube
    GRGeodesicTest,      // GR geodesic test (radial infall, circular orbits)
    GRBondi,             // GR-GSPH Bondi accretion (geodesic and sonic point)
    NSMerger2D,
    BNSCocoon1D,
    BNSCocoon2D,
    IsothermalSlab,
    PolytropicSlab,
    SinusoidalPerturbation,
    JeansInstability,
    BonnorEbertKI2000,  // K&I 2000 Bonnor-Ebert pressure-truncated sphere
    LaneEmdenKI2000,    // Lane-Emden density + K&I 2000 temperatures (NOT hydrostatic)
    IsothermalBonnorEbert,  // True isothermal Bonnor-Ebert sphere (self-gravitating equilibrium)
    HVCCIsothermal10K,  // 10K isothermal HVCC for IMBH-cloud interaction (Oka 2017)
    MHDShockTube1,      // MHD Shock Tube 1: Dai-Woodward (Iwasaki & Inutsuka 2011)
    MHDShockTube2,      // MHD Shock Tube 2: Strong shock (Iwasaki & Inutsuka 2011)
    MHDOrszagTang,      // MHD Orszag-Tang vortex: 2D turbulent MHD benchmark
    SRMHDBalsara1,      // SRMHD Balsara Test 1: Relativistic MHD shock
    UniformCloud,       // Uniform density cloud for IMBH-cloud tidal interaction
    DoNotUse,
};

class Solver {
    std::shared_ptr<SPHParameters>  m_param;
    std::shared_ptr<OutputManager>  m_output_manager;
    UnitSystem                      m_units;
    std::string                     m_output_dir;
    std::shared_ptr<Simulation>     m_sim;
    int                             m_snapshot_counter;
    
    // modules
    std::shared_ptr<Module> m_timestep;
    std::shared_ptr<Module> m_pre;
    std::shared_ptr<Module> m_fforce;
    std::shared_ptr<Module> m_gforce;
    
    // External forces (IMBH, moving potentials, etc.)
    std::shared_ptr<external_forces::PointMassBlackHole> m_external_bh;
    bool m_use_external_bh;
    
    // Cloud initial conditions (for shifting particles after snapshot load)
    vec_t m_cloud_initial_position;
    vec_t m_cloud_initial_velocity;
    bool m_has_cloud_initial_conditions;
    
    // relaxation
    std::shared_ptr<LaneEmdenRelaxation> m_lane_emden_relax;
    std::shared_ptr<PolytropicSlabRelaxation> m_polytropic_slab_relax;
    std::shared_ptr<PolytropicSlab2DRelaxation> m_polytropic_slab_2d_relax;
    std::shared_ptr<KoyamaInutsukaRelaxation> m_ki2000_relax;
    std::shared_ptr<IsothermalRelaxation> m_isothermal_relax;
    bool m_use_relaxation;
    int m_relaxation_steps;
    int m_relaxation_output_freq;  // Output frequency during relaxation
    bool m_relaxation_only;  // If true, exit after relaxation without simulation

    // GLASS pre-relaxation to uniformize particle spacing
    bool m_use_glass_relaxation;
    int m_glass_relaxation_steps;
    int m_glass_target_neighbors;
    
    // Resume configuration (JSON as SSOT - snapshot provides particle data only)
    bool m_resume_from_checkpoint;
    std::string m_checkpoint_file;
    bool m_reset_time_on_resume;  // Reset simulation time to 0 when resuming
    real m_ghost_density_factor;  // Factor to multiply ghost particle density (default: 1.0)
    real m_spherical_boundary_radius; // Reflecting spherical boundary radius (0 = disabled)

    void read_parameterfile(const char * filename);
    void make_initial_condition();
    void initialize();
    void predict();
    void correct();
    void integrate();
    void update_ghost_particles();  // Mirror ghost particle properties from real particles
    void create_periodic_ghost_particles();  // Create ghost particles for periodic boundaries
    void remove_periodic_ghost_particles();  // Remove ghost particles before next cycle
    
    // Helper method to compute total energies
    void compute_total_energies(real& kinetic, real& thermal, real& potential) const;

    // for sample
    Sample                                      m_sample;
    std::unordered_map<std::string, boost::any> m_sample_parameters;

    void make_shock_tube();
    void make_shock_tube_2d();
    void make_shock_tube_3d();
    void make_vacuum();
    void make_strong_shock();
    void make_pressure_equilibrium();
    void make_gresho_chan_vortex();
    void make_pairing_instability();
    void make_hydrostatic();
    void make_ism_cooling_1d();
    void make_khi();
    void make_evrard();
    void make_evrard_cold_collapse();
    void make_lane_emden();
    void make_lane_emden_cylinder();    // 3D cylindrical Lane-Emden
    void make_polytropic_slab_2d();     // 2D planar Lane-Emden slab
    void make_sedov();
    void make_sr_sod();
    void make_sr_tangent_velocity();
    void make_sr_rosswog();
    void make_gr_schwarzschild_shock();  // GR-GSPH Schwarzschild radial shock tube
    void make_gr_geodesic_test();        // GR geodesic test (radial infall, circular orbits)
    void make_gr_bondi();                // GR-GSPH Bondi accretion
    void make_ns_merger_2d();
    void make_bns_cocoon_1d();
    void make_bns_cocoon_2d();
    void make_isothermal_slab();
    void make_polytropic_slab();
    void make_sinusoidal_perturbation();
    void make_jeans_instability();
    void make_bonnor_ebert_ki2000();  // K&I 2000 pressure-truncated sphere
    void make_lane_emden_ki2000();    // Lane-Emden + K&I temperatures
    void make_isothermal_bonnor_ebert();  // True isothermal BE sphere (self-gravitating)
    void make_hvcc_isothermal_10k();  // 10K isothermal HVCC (Oka 2017)
    void make_mhd_shock_tube_1();     // MHD Shock Tube 1: Dai-Woodward
    void make_mhd_shock_tube_2();     // MHD Shock Tube 2: Strong shock
    void make_mhd_orszag_tang();      // MHD Orszag-Tang vortex: 2D turbulent MHD
    void make_srmhd_balsara_1();      // SRMHD Balsara Test 1: Relativistic MHD shock
    void make_uniform_cloud();        // Uniform density cloud for IMBH-cloud

public:
    Solver(int argc, char * argv[]);
    void run();
};

}