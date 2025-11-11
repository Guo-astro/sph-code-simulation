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
	@echo "  make lane_emden_help    # Lane-Emden polytrope targets"
	@echo "  make shock_tube_help    # Shock tube test targets"
	@echo ""
	@false

help: all

# Lane-Emden preset-based system
-include lane_emden/Makefile.lane_emden

# Shock Tube preset-based system
-include sample/shock_tube/Makefile.shock_tube
