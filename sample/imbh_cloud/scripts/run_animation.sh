#!/bin/bash
# Run animation script in a way that won't be interrupted

cd /Users/guo/Downloads/sphcode

echo "Starting animation generation at $(date)"
echo "This will take several minutes for 162 frames..."

python3 sample/imbh_cloud/scripts/animate_single_run.py \
    --input-dir sample/imbh_cloud/results/imbh_relaxed_2k_b3pc_nocool \
    --output-dir sample/imbh_cloud/results/imbh_relaxed_2k_b3pc_nocool \
    --mode tidal --bh-mass 100000 --fps 15 --dpi 80

echo ""
echo "Animation finished at $(date)"
ls -lh sample/imbh_cloud/results/imbh_relaxed_2k_b3pc_nocool/*.gif
