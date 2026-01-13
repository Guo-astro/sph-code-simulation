#!/bin/bash
# Compile LaTeX paper
# Requires: texlive or mactex with revtex4-2 package

cd "$(dirname "$0")"

echo "Compiling TDE_Theory_Paper.tex..."

# First pass
pdflatex -interaction=nonstopmode TDE_Theory_Paper.tex

# BibTeX (if needed)
# bibtex TDE_Theory_Paper

# Second pass for references
pdflatex -interaction=nonstopmode TDE_Theory_Paper.tex

# Third pass for final
pdflatex -interaction=nonstopmode TDE_Theory_Paper.tex

echo ""
echo "Done! Output: TDE_Theory_Paper.pdf"
echo ""

# Clean auxiliary files (optional)
# rm -f *.aux *.log *.out *.toc *.bbl *.blg
