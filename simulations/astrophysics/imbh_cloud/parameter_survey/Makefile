# IMBH-Cloud Parameter Survey Makefile
# =====================================

SURVEY_DIR := $(shell pwd)
SPH_BIN := $(SURVEY_DIR)/../../../../build/sph
PORT := 8080

.PHONY: help serve viz run-0.13pc run-0.4pc clean

help:
	@echo "IMBH-Cloud Parameter Survey"
	@echo "============================"
	@echo ""
	@echo "Visualization:"
	@echo "  make serve      - Start HTTP server for visualization (port $(PORT))"
	@echo "  make viz        - Alias for 'make serve'"
	@echo ""
	@echo "Simulations:"
	@echo "  make run-0.13pc - Run r_p=0.13pc simulations (adiabatic + cooling)"
	@echo "  make run-0.4pc  - Run r_p=0.4pc simulations (adiabatic + cooling)"
	@echo ""
	@echo "Utilities:"
	@echo "  make ic R=0.5   - Compute initial conditions for r_p=0.5pc"
	@echo "  make clean      - Remove result files"
	@echo ""
	@echo "After 'make serve', open: http://localhost:$(PORT)/viz/imbh_cloud_viz.html"

# Start HTTP server for visualization
serve:
	@echo "Starting HTTP server on port $(PORT)..."
	@echo "Open: http://localhost:$(PORT)/viz/imbh_cloud_viz.html"
	@echo "Press Ctrl+C to stop"
	@cd $(SURVEY_DIR) && python3 -m http.server $(PORT)

viz: serve

# Run r_p=0.13pc simulations
run-0.13pc:
	@echo "Running r_p=0.13pc adiabatic simulation..."
	$(SPH_BIN) parabolic_rp0.13pc/config/adiabatic.json &
	@echo "Running r_p=0.13pc cooling simulation..."
	$(SPH_BIN) parabolic_rp0.13pc/config/cooling.json &
	@echo "Simulations started in background. Check results in parabolic_rp0.13pc/*/results/"

# Run r_p=0.4pc simulations
run-0.4pc:
	@echo "Running r_p=0.4pc adiabatic simulation..."
	$(SPH_BIN) parabolic_rp0.4pc/config/adiabatic.json &
	@echo "Running r_p=0.4pc cooling simulation..."
	$(SPH_BIN) parabolic_rp0.4pc/config/cooling.json &
	@echo "Simulations started in background. Check results in parabolic_rp0.4pc/*/results/"

# Compute initial conditions for custom r_p
ic:
ifndef R
	@echo "Usage: make ic R=<pericenter_distance_in_pc>"
	@echo "Example: make ic R=0.5"
else
	python3 compute_ic.py --r_peri $(R) --start_distance 15.0
endif

# Clean result files (be careful!)
clean:
	@echo "This will remove all simulation results. Are you sure? [y/N]"
	@read ans && [ "$$ans" = "y" ] && \
		rm -rf parabolic_rp*/adiabatic/results/*.csv parabolic_rp*/cooling/results/*.csv || \
		echo "Aborted."
