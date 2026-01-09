# =============================================================================
# IMBH-Cloud Interaction Simulation (Option A: Atomic CNM)
# Three-Phase Workflow: Relaxation -> Hydrostatic Test -> IMBH Flyby
# =============================================================================
#
# WORKFLOW:
#   make phase1   -> Phase 1: Relax BE sphere with ghost envelope
#   make phase2   -> Phase 2: Test hydrostatic stability (no IMBH)
#   make phase3   -> Phase 3: IMBH tidal interaction
#   make plot     -> Generate analysis plots with analytical overlay
#
# Option A Parameters (Atomic CNM on K&I 2000 curve):
#   mu = 1.27, T = 9 K, M = 40 Msun, xi_s = 6.45 (critical)
#
# =============================================================================

ROOT_DIR := $(shell cd ../../.. && pwd)
BUILD_DIR := $(ROOT_DIR)/build_test
SPH_EXEC := $(BUILD_DIR)/sph
CONFIG_DIR := config
RESULTS_DIR := results
SCRIPTS_DIR := scripts

# Option A configs
CONFIG_PHASE1 := $(CONFIG_DIR)/presets/phase1_relaxation.json
CONFIG_PHASE2 := $(CONFIG_DIR)/presets/phase2_hydrostatic.json
CONFIG_PHASE3 := $(CONFIG_DIR)/presets/phase3_imbh_interaction.json

# Results directories
RESULTS_PHASE1 := $(RESULTS_DIR)/phase1_relaxation
RESULTS_PHASE2 := $(RESULTS_DIR)/phase2_hydrostatic
RESULTS_PHASE3 := $(RESULTS_DIR)/phase3_imbh_interaction

.PHONY: all
all: info

# =============================================================================
# Build
# =============================================================================

.PHONY: build
build:
	@echo "Building SPH code..."
	@cd $(BUILD_DIR) && make -j4 sph

# =============================================================================
# Phase 1: Relaxation with Ghost Envelope
# =============================================================================

.PHONY: phase1 relax
phase1: build
	@echo ""
	@echo "=== Phase 1: BE Sphere Relaxation (Option A) ==="
	@echo "mu = 1.27 (atomic), T = 9 K, M = 40 Msun"
	@echo "Creating relaxed hydrostatic equilibrium with ghost envelope"
	@echo ""
	@mkdir -p $(RESULTS_PHASE1)
	$(SPH_EXEC) $(CONFIG_PHASE1)
	@echo ""
	@echo "Phase 1 complete. Run: make plot-phase1"

relax: phase1

# =============================================================================
# Phase 2: Hydrostatic Stability Test
# =============================================================================

.PHONY: phase2 hydro
phase2: build
	@echo ""
	@echo "=== Phase 2: Hydrostatic Stability Test ==="
	@echo "Free evolution without IMBH to verify equilibrium"
	@echo ""
	@if [ ! -d $(RESULTS_PHASE1) ]; then \
		echo "ERROR: Run 'make phase1' first"; exit 1; \
	fi
	@mkdir -p $(RESULTS_PHASE2)
	$(SPH_EXEC) $(CONFIG_PHASE2)
	@echo ""
	@echo "Phase 2 complete. Run: make plot-phase2"

hydro: phase2

# =============================================================================
# Phase 3: IMBH Tidal Interaction
# =============================================================================

.PHONY: phase3 interact
phase3: build
	@echo ""
	@echo "=== Phase 3: IMBH Tidal Interaction ==="
	@echo "M_BH = 10^5 Msun, d = 10 pc, v = 100 km/s"
	@echo ""
	@if [ ! -d $(RESULTS_PHASE1) ]; then \
		echo "ERROR: Run 'make phase1' first"; exit 1; \
	fi
	@mkdir -p $(RESULTS_PHASE3)
	$(SPH_EXEC) $(CONFIG_PHASE3)
	@echo ""
	@echo "Phase 3 complete. Run: make plot-phase3"

interact: phase3

# =============================================================================
# Plotting (with analytical BE profile overlay)
# =============================================================================

.PHONY: plot plot-phase1 plot-phase2 plot-all
plot-phase1:
	@echo "Generating Phase 1 analysis plots..."
	@python3 $(SCRIPTS_DIR)/plot_be_sphere_analysis.py \
		--results $(RESULTS_PHASE1) \
		--config $(CONFIG_PHASE1)

plot-phase2:
	@echo "Generating Phase 2 analysis plots..."
	@python3 $(SCRIPTS_DIR)/plot_be_sphere_analysis.py \
		--results $(RESULTS_PHASE2) \
		--config $(CONFIG_PHASE2)

plot-all: plot-phase1 plot-phase2
	@echo "All plots generated."

plot: plot-phase1

# =============================================================================
# Full Workflow
# =============================================================================

.PHONY: run run-quick
run: phase1 phase2 phase3 plot-all
	@echo ""
	@echo "=== All phases complete ==="

run-quick:
	@echo "Quick test mode: reduced steps"
	@$(MAKE) phase1 RELAXATION_STEPS=1000

# =============================================================================
# Utilities
# =============================================================================

.PHONY: clean clean-phase1 clean-phase2 clean-phase3
clean-phase1:
	rm -rf $(RESULTS_PHASE1)

clean-phase2:
	rm -rf $(RESULTS_PHASE2)

clean-phase3:
	rm -rf $(RESULTS_PHASE3)

clean: clean-phase1 clean-phase2 clean-phase3
	@echo "All results cleaned."

.PHONY: info help
info:
	@echo ""
	@echo "IMBH-Cloud Simulation - Option A (Atomic CNM)"
	@echo "============================================="
	@echo ""
	@echo "Parameters: mu=1.27, T=9K, M=40 Msun, xi_s=6.45"
	@echo ""
	@echo "Phases:"
	@echo "  make phase1     Phase 1: Relax BE sphere (ghost envelope)"
	@echo "  make phase2     Phase 2: Hydrostatic stability test"
	@echo "  make phase3     Phase 3: IMBH tidal interaction"
	@echo ""
	@echo "Plotting:"
	@echo "  make plot-phase1  Density profile with analytical BE curve"
	@echo "  make plot-phase2  Hydrostatic evolution analysis"
	@echo "  make plot-all     All plots"
	@echo ""
	@echo "Workflow:"
	@echo "  make run          Run all phases + plots"
	@echo "  make clean        Delete all results"
	@echo ""

help: info
