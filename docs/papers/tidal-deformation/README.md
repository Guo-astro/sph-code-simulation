# Tidal Disruption Event Papers Collection

This directory contains key research papers on tidal disruption events (TDEs), pancaking, and shock formation relevant to the IMBH-cloud interaction simulations.

## Downloaded Papers

### ✅ PDFs Available (5 papers)

| Paper | File | Size | arXiv ID |
|-------|------|------|----------|
| Bonnerot & Lu (2022) - Nozzle Shock | `bonnerot_2022_nozzle_shock.pdf` | 2.6 MB | 2106.01376 |
| Coughlin et al. (2020) - Pancakes | `coughlin_2020_pancakes.pdf` | 2.7 MB | 2008.12304 |
| Guillochon & Ramirez-Ruiz (2013) | `guillochon_2013_tidal_disruption.pdf` | 4.7 MB | 1206.2350 |
| Stone et al. (2013) - Simulations | `stone_2013_tde_simulations.pdf` | 1.2 MB | 1210.3374 |
| Lodato et al. (2021) - Review | `lodato_2021_tde_review.pdf` | 498 KB | 2012.01817 |

### ✅ LaTeX Sources Available (5 papers)

All papers have extracted TeX source in `sources/[paper_name]_source/` directories.

## Directory Structure

```
tidal-deformation/
├── README.md                    (this file)
├── DOWNLOAD_GUIDE.md            (detailed download instructions)
├── download_all.sh              (automated download script)
├── pdfs/                        (PDF files)
│   ├── bonnerot_2022_nozzle_shock.pdf
│   ├── coughlin_2020_pancakes.pdf
│   ├── guillochon_2013_tidal_disruption.pdf
│   ├── lodato_2021_tde_review.pdf
│   └── stone_2013_tde_simulations.pdf
└── sources/                     (LaTeX source files)
    ├── bonnerot_2022_source/
    ├── coughlin_2020_source/
    ├── guillochon_2013_source/
    ├── lodato_2021_source/
    └── stone_2013_source/
```

## Key Physics Covered

### 1. **Bonnerot & Lu (2022)** - The Nozzle Shock
- Vertical compression at pericenter
- Nozzle shock formation
- Energy dissipation in TDEs
- **Relevance**: Provides the rigorous treatment of pancaking shocks

### 2. **Coughlin et al. (2020)** - Persistence of Pancakes
- Self-gravity revival post-pericenter
- In-plane caustic formation
- Pancake morphology
- **Relevance**: Explains why pancaking persists even for distant encounters

### 3. **Guillochon & Ramirez-Ruiz (2013)** - Impact Parameter Study
- Systematic study of impact parameter effects
- Stellar structure dependence
- Fallback rates
- **Relevance**: Direct analog to our cloud-IMBH parameter space study

### 4. **Stone et al. (2013)** - Hydrodynamic Simulations
- Strong compression consequences
- Numerical methods for TDEs
- Shock heating calculations
- **Relevance**: Benchmark for our SPH simulations

### 5. **Lodato et al. (2021)** - Comprehensive Review
- Modern overview of TDE physics
- Observational signatures
- Theoretical framework
- **Relevance**: Comprehensive reference covering Rees (1988) physics and beyond

## Missing Paper

### ❌ Rees (1988) - Original TDE Paper
**Why not included**: Published in Nature before arXiv existed (1991)

**Alternative access**:
- ADS: https://ui.adsabs.harvard.edu/abs/1988Natur.333..523R
- Nature (paywall): https://www.nature.com/articles/333523a0
- **Recommended**: Read Lodato et al. (2021) review which covers the same physics

## Usage

### Read a paper
```bash
# PDFs
open pdfs/bonnerot_2022_nozzle_shock.pdf

# Or use your preferred PDF reader
```

### Access LaTeX source
```bash
# Navigate to source directory
cd sources/bonnerot_2022_source

# Main TeX file is usually named ms.tex, main.tex, or arxiv.tex
ls *.tex
```

### Re-download all papers
```bash
./download_all.sh
```

## Citation Information

See `DOWNLOAD_GUIDE.md` for complete BibTeX entries for all papers.

## Relevance to IMBH-Cloud Simulations

These papers establish the theoretical framework for:

1. **Shock velocity formula**: v_shock = R_cloud × √(G M_BH / b³)
2. **Pancaking physics**: Vertical compression perpendicular to orbital plane
3. **Observable signatures**: When molecules survive vs dissociate
4. **Numerical methods**: SPH techniques for tidal disruption

All physics is detailed in:
```
/Users/guo/Downloads/sphcode/sample/imbh_cloud/docs/RESEARCH_SETUP.md
```
Section A.9: "Shock Heating and Molecular Dissociation"

---

**Last updated**: 2025-12-02
**Total size**: ~25 MB (PDFs + sources)
