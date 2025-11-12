# SPH Code Simulation - Main Makefile
#
# Build system: Use CMake in build/ directory
#   cd build && cmake .. && make -j8
#
# Simulation targets: Preset-based workflows for different samples

.PHONY: all help

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
	@echo ""
	@false

help: all

# Lane-Emden preset-based system
-include lane_emden/Makefile.lane_emden

# Shock Tube preset-based system
-include sample/shock_tube/Makefile.shock_tube

# Shock Tube preset-based system
-include sample/shock_tube_2d/Makefile.shock_tube_2d

# Pairing Instability preset-based system
-include sample/pairing_instability/Makefile.pairing_instability

# Hydrostatic Test preset-based system
-include sample/hydrostatic/Makefile.hydrostatic

# Kelvin-Helmholtz Instability preset-based system
-include sample/khi/Makefile.khi

# Gresho-Chan Vortex preset-based system
-include sample/gresho_chan_vortex/Makefile.gresho

# Sedov Blast Wave preset-based system
-include sample/sedov/Makefile.sedov

# Vacuum Test preset-based system
-include sample/vacuum/Makefile.vacuum

# Strong Shock Test preset-based system
-include sample/strong_shock/Makefile.strong_shock
