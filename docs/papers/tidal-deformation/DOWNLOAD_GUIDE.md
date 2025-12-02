# Tidal Disruption Papers - Download Guide

This directory contains papers referenced in the IMBH-cloud interaction research documentation.

## Papers with arXiv Access (PDF + Source Available)

### 1. Bonnerot & Lu (2022) - The Nozzle Shock in TDEs
**arXiv ID**: 2106.01376
**Journal**: MNRAS 511, 2147 (2022)

**Download commands:**
```bash
# PDF
wget https://arxiv.org/pdf/2106.01376.pdf -O bonnerot_2022_nozzle_shock.pdf

# Source (tar.gz with TeX files)
wget https://arxiv.org/e-print/2106.01376 -O bonnerot_2022_source.tar.gz
tar -xzf bonnerot_2022_source.tar.gz
```

**Links:**
- arXiv: https://arxiv.org/abs/2106.01376
- Journal: https://academic.oup.com/mnras/article/511/2/2147/6516431

---

### 2. Coughlin et al. (2020) - Persistence of Pancakes
**arXiv ID**: 2008.12304
**Journal**: ApJ Letters 900, L45 (2020)

**Download commands:**
```bash
# PDF
wget https://arxiv.org/pdf/2008.12304.pdf -O coughlin_2020_pancakes.pdf

# Source
wget https://arxiv.org/e-print/2008.12304 -O coughlin_2020_source.tar.gz
tar -xzf coughlin_2020_source.tar.gz
```

**Links:**
- arXiv: https://arxiv.org/abs/2008.12304
- Journal: https://iopscience.iop.org/article/10.3847/2041-8213/abb2ad

---

### 3. Guillochon & Ramirez-Ruiz (2013) - Impact Parameter & Stellar Structure
**arXiv ID**: 1206.2350
**Journal**: ApJ 767, 25 (2013)

**Download commands:**
```bash
# PDF
wget https://arxiv.org/pdf/1206.2350.pdf -O guillochon_2013_tidal_disruption.pdf

# Source
wget https://arxiv.org/e-print/1206.2350 -O guillochon_2013_source.tar.gz
tar -xzf guillochon_2013_source.tar.gz
```

**Links:**
- arXiv: https://arxiv.org/abs/1206.2350
- DOI: 10.1088/0004-637X/767/1/25

---

## Papers NOT on arXiv (Journal/Publisher Access Only)

### 4. Rees (1988) - Seminal TDE Paper
**Journal**: Nature 333, 523-528 (1988)
**DOI**: 10.1038/333523a0

**Note**: This paper predates arXiv (founded 1991). Access requires:
- Institutional subscription to Nature
- Purchase from publisher
- Request from author via ResearchGate

**Links:**
- ADS: https://ui.adsabs.harvard.edu/abs/1988Natur.333..523R/abstract
- Nature: https://www.nature.com/articles/333523a0

**Alternative**: Many textbooks and review articles summarize the key physics from this paper.

---

## Additional Relevant Papers (Recommended)

### 5. Stone et al. (2013) - Hydrodynamic Simulations of TDEs
**arXiv ID**: 1210.3374
**Journal**: MNRAS 435, 1809 (2013)

```bash
wget https://arxiv.org/pdf/1210.3374.pdf -O stone_2013_tde_simulations.pdf
wget https://arxiv.org/e-print/1210.3374 -O stone_2013_source.tar.gz
```

---

### 6. Lodato et al. (2021) - Review: The Process of Stellar Tidal Disruption
**arXiv ID**: 2012.01817
**Journal**: Space Science Reviews 217, 83 (2021)

```bash
wget https://arxiv.org/pdf/2012.01817.pdf -O lodato_2021_tde_review.pdf
wget https://arxiv.org/e-print/2012.01817 -O lodato_2021_source.tar.gz
```

**Links:**
- arXiv: https://arxiv.org/abs/2012.01817
- Journal: https://link.springer.com/article/10.1007/s11214-021-00818-7

---

## Quick Download All Available Papers

Run this script to download all papers with arXiv access:

```bash
#!/bin/bash
cd /Users/guo/Downloads/sphcode/docs/papers/tidal-deformation

# Create subdirectories
mkdir -p pdfs sources

# Download PDFs
wget https://arxiv.org/pdf/2106.01376.pdf -O pdfs/bonnerot_2022_nozzle_shock.pdf
wget https://arxiv.org/pdf/2008.12304.pdf -O pdfs/coughlin_2020_pancakes.pdf
wget https://arxiv.org/pdf/1206.2350.pdf -O pdfs/guillochon_2013_tidal_disruption.pdf
wget https://arxiv.org/pdf/1210.3374.pdf -O pdfs/stone_2013_tde_simulations.pdf
wget https://arxiv.org/pdf/2012.01817.pdf -O pdfs/lodato_2021_tde_review.pdf

# Download sources
wget https://arxiv.org/e-print/2106.01376 -O sources/bonnerot_2022_source.tar.gz
wget https://arxiv.org/e-print/2008.12304 -O sources/coughlin_2020_source.tar.gz
wget https://arxiv.org/e-print/1206.2350 -O sources/guillochon_2013_source.tar.gz
wget https://arxiv.org/e-print/1210.3374 -O sources/stone_2013_source.tar.gz
wget https://arxiv.org/e-print/2012.01817 -O sources/lodato_2021_source.tar.gz

echo "Download complete! PDFs in pdfs/, sources in sources/"
```

Save this as `download_all.sh` and run:
```bash
chmod +x download_all.sh
./download_all.sh
```

---

## Citation Information

### BibTeX Entries

```bibtex
@article{Rees1988,
  author = {Rees, Martin J.},
  title = {Tidal disruption of stars by black holes of $10^6$--$10^8$ solar masses in nearby galaxies},
  journal = {Nature},
  volume = {333},
  pages = {523--528},
  year = {1988},
  doi = {10.1038/333523a0}
}

@article{Bonnerot2022,
  author = {Bonnerot, Cl{\'e}ment and Lu, Wenbin},
  title = {The nozzle shock in tidal disruption events},
  journal = {MNRAS},
  volume = {511},
  pages = {2147--2169},
  year = {2022},
  doi = {10.1093/mnras/stac146},
  archivePrefix = {arXiv},
  eprint = {2106.01376}
}

@article{Coughlin2020,
  author = {Coughlin, Eric R. and Nixon, C. J. and Begelman, Mitchell C. and Armitage, Philip J.},
  title = {The Persistence of Pancakes and the Revival of Self-gravity in Tidal Disruption Events},
  journal = {ApJ Letters},
  volume = {900},
  pages = {L45},
  year = {2020},
  doi = {10.3847/2041-8213/abb2ad},
  archivePrefix = {arXiv},
  eprint = {2008.12304}
}

@article{Guillochon2013,
  author = {Guillochon, James and Ramirez-Ruiz, Enrico},
  title = {Hydrodynamical Simulations to Determine the Feeding Rate of Black Holes by the Tidal Disruption of Stars: The Importance of the Impact Parameter and Stellar Structure},
  journal = {ApJ},
  volume = {767},
  pages = {25},
  year = {2013},
  doi = {10.1088/0004-637X/767/1/25},
  archivePrefix = {arXiv},
  eprint = {1206.2350}
}

@article{Lodato2021,
  author = {Lodato, G. and Rossi, E. M. and others},
  title = {The Process of Stellar Tidal Disruption by Supermassive Black Holes},
  journal = {Space Science Reviews},
  volume = {217},
  pages = {83},
  year = {2021},
  doi = {10.1007/s11214-021-00818-7},
  archivePrefix = {arXiv},
  eprint = {2012.01817}
}
```

---

## Summary

- ✅ **5 papers** available on arXiv (PDF + TeX source)
- ❌ **1 paper** (Rees 1988) requires institutional access
- 📚 **Total**: 6 key papers on tidal disruption physics

For the Rees (1988) paper, consider using the review by Lodato et al. (2021) which comprehensively covers the same physics with modern context.
