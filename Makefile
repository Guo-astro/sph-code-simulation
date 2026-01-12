# SPH Code Simulation - Main Makefile
#
# Build system: Use CMake in build/ directory
#   cd build && cmake .. && make -j8
#
# Simulation targets: Preset-based workflows for different samples

.PHONY: all help viz viz_export viz_server

# Default target
all:
	@echo "=========================================="
	@echo "SPH Code Simulation"
	@echo "=========================================="
	@echo ""
	@echo "Build with CMake:"
	@echo "  cd build && make -j8"
	@echo ""
	@echo "Available simulation workflows:"
	@echo "  make lane_emden_help          # Lane-Emden polytrope targets"
	@echo "  make shock_tube_help          # Shock tube test targets"
	@echo "  make pairing_help             # Pairing instability targets"
	@echo "  make hydrostatic_help         # Hydrostatic test targets"
	@echo "  make khi_help                 # Kelvin-Helmholtz instability targets"
	@echo "  make gresho_help              # Gresho-Chan vortex targets"
	@echo "  make sedov_help               # Sedov blast wave targets"
	@echo "  make vacuum_help              # Vacuum test targets"
	@echo "  make strong_shock_help        # Strong shock test targets"
	@echo "  make sr_help                  # SR-GSPH relativistic test targets"
	@echo "  make imbh_help                # IMBH-cloud tidal disruption targets"
	@echo ""
	@echo "Grad-h Study Test Suites:"
	@echo "  make gradh_2d_help            # 2D Planar Lane-Emden grad-h study"
	@echo "  make gradh_3d_help            # 3D Cylinder Lane-Emden grad-h study"
	@echo ""
	@echo "Visualization:"
	@echo "  make viz SIM=<path>           # Export data and start viz server"
	@echo "  make viz_export SIM=<path>    # Export simulation data only"
	@echo "  make viz_server               # Start visualization server only"
	@echo ""
	@echo "Quick IMBH shortcuts:"
	@echo "  make imbh_relax_2k            # Testing: 2k particle relaxation"
	@echo "  make imbh_relax_200k          # Production: 200k particle relaxation"
	@echo ""
	@false

help: all

#==============================================================================
# SPH Visualization Tool
#==============================================================================
# Usage:
#   make viz SIM=simulations/benchmarks/sedov/results/gsph_wendland
#   make viz SIM=simulations/astrophysics/imbh_cloud/results/Mc1e3_Mbh1e5_b3_v10/adiabatic_61k_gsph
#   make viz_export SIM=simulations/astrophysics/lane_emden/results/n3_gsph
#   make viz_server
#
# Options:
#   SIM         Path to simulation results directory (required for viz/viz_export)
#   STRIDE      Export every Nth frame (default: 1)
#   MAX_FRAMES  Maximum frames to export (default: all)
#==============================================================================

VIZ_DIR := tools/sph-viz
VIZ_EXPORT := python $(VIZ_DIR)/scripts/export_viz_data.py
STRIDE ?= 1
MAX_FRAMES ?=

# Export simulation data for visualization
viz_export:
ifndef SIM
	@echo "❌ Error: SIM path required"
	@echo "Usage: make viz_export SIM=<simulation_path>"
	@echo "Example: make viz_export SIM=simulations/benchmarks/sedov/results/gsph_wendland"
	@exit 1
endif
	@echo "========================================"
	@echo "📊 Exporting visualization data"
	@echo "========================================"
	@echo "Simulation: $(SIM)"
	$(VIZ_EXPORT) $(SIM) --stride $(STRIDE) $(if $(MAX_FRAMES),--max-frames $(MAX_FRAMES),)
	@echo "✓ Export complete!"

# Start visualization server
viz_server:
	@echo "========================================"
	@echo "🚀 Starting SPH Visualization Server"
	@echo "========================================"
	@echo ""
	@echo "Opening: http://localhost:3000/viz"
	@echo "Press Ctrl+C to stop the server"
	@echo ""
	@cd $(VIZ_DIR) && npm run dev

# One-shot: export data and start server
viz: viz_export viz_server

# Lane-Emden preset-based system
-include simulations/astrophysics/lane_emden/Makefile.lane_emden
-include simulations/astrophysics/lane_emden/Makefile.2d
-include simulations/astrophysics/lane_emden/Makefile.3d

# Shock Tube preset-based system
-include simulations/benchmarks/shock_tube/Makefile.shock_tube

# Shock Tube preset-based system
-include simulations/benchmarks/shock_tube_2d/Makefile.shock_tube_2d

# Pairing Instability preset-based system
-include simulations/stability/pairing_instability/Makefile.pairing_instability

# Hydrostatic Test preset-based system
-include simulations/stability/hydrostatic/Makefile.hydrostatic

# Kelvin-Helmholtz Instability preset-based system
-include simulations/stability/khi/Makefile.khi

# Gresho-Chan Vortex preset-based system
-include simulations/stability/gresho_chan_vortex/Makefile.gresho

# Sedov Blast Wave preset-based system
-include simulations/benchmarks/sedov/Makefile.sedov

# Vacuum Test preset-based system
-include simulations/benchmarks/vacuum/Makefile.vacuum

# Strong Shock Test preset-based system
-include simulations/benchmarks/strong_shock/Makefile.strong_shock

# SR-GSPH Relativistic Test preset-based system
-include simulations/relativistic/sr_sod/Makefile.sr_sod

# IMBH-Cloud Tidal Disruption preset-based system
-include simulations/astrophysics/imbh_cloud/Makefile.imbh_cloud

# Grad-h Study Test Suites (GSPH/SSPH × grad-h comparison)
-include simulations/stability/gradh_study_2d_planar/Makefile.gradh_2d_planar
-include simulations/stability/gradh_study_3d_cylinder/Makefile.gradh_3d_cylinder
