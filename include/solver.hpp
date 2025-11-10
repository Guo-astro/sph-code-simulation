#pragma once

#include <memory>
#include <unordered_map>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/any.hpp>

#include "defines.hpp"

namespace sph
{

struct SPHParameters;
class Simulation;
class Output;
class LaneEmdenRelaxation;

class Module;

enum struct Sample {
    ShockTube,
    GreshoChanVortex,
    PairingInstability,
    HydroStatic,
    KHI,
    Evrard,
    EvrardColdCollapse,
    LaneEmden,
    DoNotUse,
};

class Solver {
    std::shared_ptr<SPHParameters>  m_param;
    std::shared_ptr<Output>         m_output;
    std::string                     m_output_dir;
    std::shared_ptr<Simulation>     m_sim;

    // modules
    std::shared_ptr<Module> m_timestep;
    std::shared_ptr<Module> m_pre;
    std::shared_ptr<Module> m_fforce;
    std::shared_ptr<Module> m_gforce;
    
    // relaxation
    std::shared_ptr<LaneEmdenRelaxation> m_lane_emden_relax;
    bool m_use_relaxation;
    int m_relaxation_steps;
    int m_relaxation_output_freq;  // Output frequency during relaxation
    bool m_relaxation_only;  // If true, exit after relaxation without simulation

    void read_parameterfile(const char * filename);
    void make_initial_condition();
    void initialize();
    void predict();
    void correct();
    void integrate();

    // for sample
    Sample                                      m_sample;
    std::unordered_map<std::string, boost::any> m_sample_parameters;

    void make_shock_tube();
    void make_gresho_chan_vortex();
    void make_pairing_instability();
    void make_hydrostatic();
    void make_khi();
    void make_evrard();
    void make_evrard_cold_collapse();
    void make_lane_emden();

public:
    Solver(int argc, char * argv[]);
    void run();
};

}