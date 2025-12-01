#!/bin/bash
# Run script for Figure 1 reproduction
# Koyama & Inutsuka (2000) ISM chemistry network

set -e

echo "========================================"
echo "Figure 1 Reproduction"
echo "Koyama & Inutsuka (2000, ApJ 533, 793)"
echo "========================================"
echo ""

# Check Python
if ! command -v python3 &> /dev/null; then
    echo "❌ Error: python3 not found"
    exit 1
fi

echo "✓ Python found: $(python3 --version)"
echo ""

# Check dependencies
echo "Checking dependencies..."
python3 -c "import numpy; print('✓ numpy:', numpy.__version__)" || {
    echo "❌ numpy not found. Install with: pip install numpy"
    exit 1
}

python3 -c "import scipy; print('✓ scipy:', scipy.__version__)" || {
    echo "❌ scipy not found. Install with: pip install scipy"
    exit 1
}

python3 -c "import matplotlib; print('✓ matplotlib:', matplotlib.__version__)" || {
    echo "❌ matplotlib not found. Install with: pip install matplotlib"
    exit 1
}

echo ""
echo "All dependencies satisfied!"
echo ""

# Create results directory
mkdir -p ../results

echo "Running equilibrium solver and plotting..."
echo ""

# Run the main script
python3 reproduce_figure1.py

echo ""
echo "========================================"
echo "✓ Figure 1 reproduction complete!"
echo "========================================"
echo ""
echo "Output files:"
echo "  - ../results/figure1_reproduction.png (4-panel combined)"
echo "  - ../results/f1a_reproduction.png (Temperature & Pressure)"
echo "  - ../results/f1b_reproduction.png (Chemical fractions)"
echo "  - ../results/f1c_reproduction.png (Heating & Cooling rates)"
echo "  - ../results/f1d_reproduction.png (Timescales)"
echo ""
echo "Compare with original figures in:"
echo "  - ../../../docs/papers/cooling-heating/f1{a,b,c,d}.ps"
echo ""
