#!/bin/bash
# Download all available tidal disruption papers from arXiv
# Usage: chmod +x download_all.sh && ./download_all.sh

set -e  # Exit on error

BASEDIR="/Users/guo/Downloads/sphcode/docs/papers/tidal-deformation"
cd "$BASEDIR"

echo "================================================"
echo "Downloading Tidal Disruption Papers from arXiv"
echo "================================================"
echo ""

# Create subdirectories
mkdir -p pdfs sources

echo "[1/5] Downloading Bonnerot & Lu (2022) - Nozzle Shock..."
curl -L https://arxiv.org/pdf/2106.01376.pdf -o pdfs/bonnerot_2022_nozzle_shock.pdf
curl -L https://arxiv.org/e-print/2106.01376 -o sources/bonnerot_2022_source.tar.gz

echo "[2/5] Downloading Coughlin et al. (2020) - Pancakes & Self-Gravity..."
curl -L https://arxiv.org/pdf/2008.12304.pdf -o pdfs/coughlin_2020_pancakes.pdf
curl -L https://arxiv.org/e-print/2008.12304 -o sources/coughlin_2020_source.tar.gz

echo "[3/5] Downloading Guillochon & Ramirez-Ruiz (2013) - Impact Parameter..."
curl -L https://arxiv.org/pdf/1206.2350.pdf -o pdfs/guillochon_2013_tidal_disruption.pdf
curl -L https://arxiv.org/e-print/1206.2350 -o sources/guillochon_2013_source.tar.gz

echo "[4/5] Downloading Stone et al. (2013) - Hydrodynamic Simulations..."
curl -L https://arxiv.org/pdf/1210.3374.pdf -o pdfs/stone_2013_tde_simulations.pdf
curl -L https://arxiv.org/e-print/1210.3374 -o sources/stone_2013_source.tar.gz

echo "[5/5] Downloading Lodato et al. (2021) - Comprehensive Review..."
curl -L https://arxiv.org/pdf/2012.01817.pdf -o pdfs/lodato_2021_tde_review.pdf
curl -L https://arxiv.org/e-print/2012.01817 -o sources/lodato_2021_source.tar.gz

echo ""
echo "================================================"
echo "✓ Download complete!"
echo "================================================"
echo ""
echo "PDFs saved to:    $BASEDIR/pdfs/"
echo "Sources saved to: $BASEDIR/sources/"
echo ""
echo "Extracting LaTeX sources..."
cd sources
for tarfile in *.tar.gz; do
    basename="${tarfile%.tar.gz}"
    mkdir -p "$basename"
    tar -xzf "$tarfile" -C "$basename" 2>/dev/null || echo "  Note: $tarfile may be a single .tex file"
done
echo ""
echo "✓ Sources extracted to subdirectories in sources/"
echo ""
echo "Summary:"
echo "  - 5 PDFs downloaded"
echo "  - 5 TeX source packages downloaded"
echo ""
echo "Note: Rees (1988) is NOT available on arXiv (predates arXiv)."
echo "      See DOWNLOAD_GUIDE.md for alternative access methods."
echo ""
