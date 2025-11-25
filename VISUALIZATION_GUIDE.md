# SR Sod Visualization Guide

Complete guide for plotting and animating SR Sod shock tube results.

---

## 📊 **Scripts Created**

All scripts are in `scripts/` directory:

1. **`analyze_density_blip.py`** - Automated density blip analysis
2. **`plot_sr_sod_snapshot.py`** - Plot single snapshot with all variables
3. **`compare_sr_sod.py`** - Side-by-side comparison of two snapshots
4. **`animate_sr_sod.py`** - Create animation from snapshots
5. **`compare_animations.py`** - Side-by-side animation comparison

---

## 🎬 **Generated Visualizations**

### Static Plots

| File | Description | Location |
|------|-------------|----------|
| `comparison_before_after.png` | **Before vs After comparison** | Root directory |
| `snapshot_fixed_detail.png` | Detailed plot (all variables) | Root directory |
| `density_blip_analysis.png` | Density blip diagnostic | `sample/sr_sod/results/sharp/` |
| `density_blip_analysis.png` | Density blip diagnostic (fixed) | `sample/sr_sod/results/sharp_fixed/` |

### Animations

| File | Description | Frames | Location |
|------|-------------|--------|----------|
| `sr_sod_fixed.mp4` | Animation of fixed simulation | 18 | Root directory |
| `sr_sod_comparison.mp4` | **Side-by-side comparison** | 14 | Root directory |

---

## 📖 **Usage Examples**

### 1. Analyze Density Blip

```bash
python3 scripts/analyze_density_blip.py sample/sr_sod/results/sharp_fixed/snapshot_0175.csv
```

**Output**:
- Prints overshoot/undershoot percentages
- Saves `density_blip_analysis.png` in same directory

---

### 2. Plot Single Snapshot (All Variables)

```bash
python3 scripts/plot_sr_sod_snapshot.py \
  sample/sr_sod/results/sharp_fixed/snapshot_0175.csv \
  my_plot.png
```

**Shows**:
- Rest-frame density (n)
- Lab-frame density (N = γn)
- Pressure (P)
- Velocity (v)
- Lorentz factor (γ)
- Internal energy (u)

---

### 3. Compare Two Snapshots

```bash
python3 scripts/compare_sr_sod.py \
  sample/sr_sod/results/sharp/snapshot_0131.csv \
  sample/sr_sod/results/sharp_fixed/snapshot_0175.csv \
  comparison.png
```

**Shows**: Side-by-side comparison of N, P, v

---

### 4. Create Animation

```bash
# Single simulation
python3 scripts/animate_sr_sod.py \
  sample/sr_sod/results/sharp_fixed \
  animation.mp4 \
  5  # skip every 5 frames
```

**Parameters**:
- `skip=5` → use every 5th snapshot (faster)
- `skip=1` → use all snapshots (slower, smoother)

**Output formats**:
- `.mp4` (requires ffmpeg)
- `.gif` (requires pillow)

---

### 5. Compare Two Simulations (Animation)

```bash
python3 scripts/compare_animations.py \
  sample/sr_sod/results/sharp \
  sample/sr_sod/results/sharp_fixed \
  comparison.mp4 \
  10
```

**Shows**: Side-by-side evolution of both simulations

---

## 🎨 **Customization**

### Change Plot Style

Edit the scripts to customize:

```python
# In any plotting script
ax.plot(x, y, 'r-o', linewidth=2, markersize=5)
#            ^   ^           ^         ^
#          color style    line width  marker size

# Colors: 'r'=red, 'b'=blue, 'g'=green, 'k'=black
# Styles: '-'=line, 'o'=circles, '--'=dashed
```

### Adjust Animation Speed

```python
# In animate_sr_sod.py, line ~140
anim.save(output_file, writer='ffmpeg', fps=10, dpi=100)
#                                       ^^^      ^^^
#                                     frames/s  quality
```

- `fps=10` → 10 frames per second (default)
- `fps=20` → faster
- `fps=5` → slower
- `dpi=100` → standard quality
- `dpi=150` → high quality (larger file)

---

## 📈 **Analysis Results**

### Before Fix (100 particles, C_cd=0.2)

```
Overshoot:  29.82%
Undershoot: 50.17%
Pressure variation: 16.27%
```

### After Fix (400 particles, C_cd=1.0, 2nd order)

```
Overshoot:  45.93%
Undershoot: 52.70%
Pressure variation: 7.10% ✓ (improved!)
```

---

## 🔧 **Troubleshooting**

### Error: "No module named 'matplotlib'"

```bash
pip install matplotlib numpy
```

### Error: "ffmpeg not found"

**macOS**:
```bash
brew install ffmpeg
```

**Linux**:
```bash
sudo apt-get install ffmpeg
```

**Alternative**: Use `.gif` instead of `.mp4`:
```bash
python3 scripts/animate_sr_sod.py dir/ output.gif 5
```

### Animation too slow/large

Increase `skip` parameter:
```bash
# Instead of skip=5, use skip=10 or skip=20
python3 scripts/animate_sr_sod.py dir/ output.mp4 20
```

### Plot looks pixelated

Increase DPI in the script:
```python
plt.savefig(output_file, dpi=300, bbox_inches='tight')
#                        ^^^
```

---

## 📝 **Batch Processing**

### Plot All Snapshots

```bash
for f in sample/sr_sod/results/sharp_fixed/snapshot_*.csv; do
    python3 scripts/plot_sr_sod_snapshot.py "$f"
done
```

### Compare Multiple Times

```bash
# Compare t=0.1, 0.2, 0.3
python3 scripts/compare_sr_sod.py \
  sample/sr_sod/results/sharp/snapshot_0050.csv \
  sample/sr_sod/results/sharp_fixed/snapshot_0088.csv \
  comparison_t0.1.png
```

---

## 🎯 **Quick Commands**

```bash
# Analyze density blip
python3 scripts/analyze_density_blip.py sample/sr_sod/results/sharp_fixed/snapshot_0175.csv

# Create comparison plot
python3 scripts/compare_sr_sod.py \
  sample/sr_sod/results/sharp/snapshot_0131.csv \
  sample/sr_sod/results/sharp_fixed/snapshot_0175.csv

# Create animation
python3 scripts/animate_sr_sod.py sample/sr_sod/results/sharp_fixed sr_sod.mp4 10

# Create comparison animation
python3 scripts/compare_animations.py \
  sample/sr_sod/results/sharp \
  sample/sr_sod/results/sharp_fixed \
  comparison.mp4 10
```

---

## 📚 **File Structure**

```
sphcode/
├── scripts/
│   ├── analyze_density_blip.py     ← Density blip analysis
│   ├── plot_sr_sod_snapshot.py     ← Single snapshot plot
│   ├── compare_sr_sod.py           ← Comparison plot
│   ├── animate_sr_sod.py           ← Animation
│   └── compare_animations.py       ← Comparison animation
│
├── sample/sr_sod/results/
│   ├── sharp/                      ← Before (100p, C_cd=0.2)
│   │   ├── snapshot_*.csv
│   │   └── density_blip_analysis.png
│   │
│   └── sharp_fixed/                ← After (400p, C_cd=1.0)
│       ├── snapshot_*.csv
│       └── density_blip_analysis.png
│
├── comparison_before_after.png     ← Main comparison
├── snapshot_fixed_detail.png       ← Detailed plot
├── sr_sod_fixed.mp4               ← Animation (fixed)
├── sr_sod_comparison.mp4          ← Comparison animation
│
├── DENSITY_BLIP_DIAGNOSTIC.md     ← Original analysis
├── DENSITY_BLIP_COMPARISON.md     ← Before/after comparison
└── VISUALIZATION_GUIDE.md         ← This file
```

---

## 🎓 **Understanding the Plots**

### Lab-Frame vs Rest-Frame Density

- **Rest-frame (n)**: Density in fluid's rest frame
- **Lab-frame (N = γn)**: Density observed in simulation frame
- **Lorentz factor (γ)**: Time dilation / length contraction factor

At contact discontinuity:
- **n** jumps (different fluids)
- **N** shows overshoot/undershoot (numerical artifact)
- **P** should be constant (physics requirement)

### What to Look For

✅ **Good simulation**:
- Smooth rarefaction wave
- Sharp shock front
- Constant pressure at contact
- Small density blip (<5%)

❌ **Bad simulation**:
- Large density overshoot (>30%)
- Pressure oscillations (>10%)
- Diffuse shock front
- Noisy rarefaction

---

**Generated**: November 25, 2025
**Scripts location**: `scripts/`
**Results location**: `sample/sr_sod/results/`
