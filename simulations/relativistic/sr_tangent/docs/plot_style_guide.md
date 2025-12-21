# Kitajima-Style Resolution Comparison Plot: Style Guide

## Overview

This document describes the detailed style specifications for the Kitajima et al. (2025) Figure 10-style resolution comparison plot. Use this guide to reproduce similar plots for other resolution studies.

---

## 1. Layout Structure

### 1.1 Grid Configuration

```
┌─────────────────────────────────────────────────────────────────┐
│                    Legend (3 columns, centered)                  │
├─────────────────────────────────────────────────────────────────┤
│              Title: Resolution Study: v^t_L = v^t_R = 0.9       │
├──────────────────┬──────────────────┬──────────────────────────┤
│ Row 0: 800+800   │                  │                          │
│ Pressure/Density │  Normal Velocity │  Tangent Velocity        │
├──────────────────┼──────────────────┼──────────────────────────┤
│ Row 1: 1600+1600 │                  │                          │
│ Pressure/Density │  Normal Velocity │  Tangent Velocity        │
├──────────────────┼──────────────────┼──────────────────────────┤
│ Row 2: 3200+1600 │                  │                          │
│ Pressure/Density │  Normal Velocity │  Tangent Velocity        │
├──────────────────┼──────────────────┼──────────────────────────┤
│ Row 3: 12800+1600│                  │                          │
│ Pressure/Density │  Normal Velocity │  Tangent Velocity        │
└──────────────────┴──────────────────┴──────────────────────────┘
```

- **Rows**: 4 (one per resolution)
- **Columns**: 3 (Pressure/Density, Normal Velocity, Tangent Velocity)
- **Figure size**: `figsize=(12, 14)`

### 1.2 Row Labels

Row labels appear as Y-axis labels on the leftmost column:
- Format: `"N_left + N_right"` with comma separators
- Examples: `"800 + 800"`, `"1,600 + 1,600"`, `"12,800 + 1,600"`
- Font: `fontsize=11, fontweight='bold'`

### 1.3 Column Titles

Column titles appear only on the first row (row 0):
- Column 0: `"Pressure & Density"`
- Column 1: `"Normal Velocity"`
- Column 2: `"Tangent Velocity"`
- Font: `fontsize=11`

### 1.4 X-axis Labels

X-axis labels appear only on the last row (row 3):
- Label: `"x"`
- Font: `fontsize=11`

---

## 2. Data Layers (Plot Order)

Data is plotted in this order (back to front):
1. **Analytical solution** (solid black lines) - background
2. **Exact Riemann solver** (x markers) - middle
3. **HLLC solver** (o markers) - foreground

---

## 3. Color Scheme

### 3.1 High-Contrast Color Palette

**Principle**: Use dark colors for Exact solver, bright colors for HLLC solver.

| Data | Exact Riemann | HLLC |
|------|---------------|------|
| Pressure | `'darkgreen'` | `'lime'` |
| Density | `'darkred'` | `'orange'` |
| v^x | `'darkblue'` | `'deepskyblue'` |
| v^t | `'purple'` | `'magenta'` |

**Analytical solution**: `'black'`

### 3.2 Colors NOT to Use

- `'cyan'` - low contrast
- Light grays - hard to see
- Similar hues for different datasets

---

## 4. Marker Styles

### 4.1 Marker Types

| Dataset | Marker | Description |
|---------|--------|-------------|
| Analytical | `-` | Solid line |
| Exact Riemann | `'x'` | Cross marker |
| HLLC | `'o'` | Circle marker (hollow) |

### 4.2 Marker Parameters

**Exact Riemann (x markers):**
```python
ax.plot(x, y, 'x', color=COLOR, ms=3, mew=0.8, alpha=0.9)
```
- `ms=3` (marker size)
- `mew=0.8` (marker edge width)
- `alpha=0.9` (transparency)

**HLLC (o markers):**
```python
ax.plot(x, y, 'o', color=COLOR, ms=4, mfc='none', mew=0.8, alpha=0.7)
```
- `ms=4` (slightly larger)
- `mfc='none'` (hollow/unfilled)
- `mew=0.8` (edge width)
- `alpha=0.7` (more transparent)

**Analytical (lines):**
```python
ax.plot(x, y, '-', color='black', lw=1.5, alpha=0.8)
```
- `lw=1.5` (line width)
- `alpha=0.8`

---

## 5. Axis Configuration

### 5.1 X-axis

- **Range**: `ax.set_xlim(-0.5, 0.5)`
- Common for all panels

### 5.2 Y-axis Ranges

| Column | Y-limits |
|--------|----------|
| Pressure/Density | `(-0.05, 1.2)` |
| Normal Velocity | `(-0.02, 0.42)` |
| Tangent Velocity | `(0.4, 1.02)` |

### 5.3 Scaling for Column 0

In the Pressure/Density column, data is scaled:
- Pressure: `P / 1000`
- Density: `n / 5`

This brings both quantities to similar ranges for visual comparison.

---

## 6. Text Annotations

### 6.1 Variable Labels

Each panel has a text label identifying the variable:

**Column 0 (Pressure/Density):**
```python
ax.text(-0.45, 1.05, 'P/1000', color='darkgreen', fontsize=10, fontweight='bold')
ax.text(-0.45, 0.25, 'n/5', color='darkred', fontsize=10, fontweight='bold')
```

**Column 1 (Normal Velocity):**
```python
ax.text(-0.45, 0.36, r'$v^x$', color='darkblue', fontsize=12, fontweight='bold')
```

**Column 2 (Tangent Velocity):**
```python
ax.text(-0.45, 0.95, r'$v^t$', color='purple', fontsize=12, fontweight='bold')
```

### 6.2 Position

- X position: `-0.45` (near left edge)
- Y position: near top of expected data range

---

## 7. Legend

### 7.1 Legend Elements

```python
from matplotlib.lines import Line2D

legend_elements = [
    Line2D([0], [0], color='black', linestyle='-', linewidth=2,
           label='Analytical'),
    Line2D([0], [0], marker='x', color='darkblue', linestyle='None',
           markersize=8, markeredgewidth=1.5, label='Exact Riemann'),
    Line2D([0], [0], marker='o', color='deepskyblue', linestyle='None',
           markersize=8, markerfacecolor='none', markeredgewidth=1.5,
           label='HLLC'),
]
```

### 7.2 Legend Placement

```python
fig.legend(handles=legend_elements, loc='upper center', ncol=3,
           fontsize=11, bbox_to_anchor=(0.5, 0.98))
```

- Position: Upper center of figure
- Columns: 3 (horizontal layout)
- Anchor: `(0.5, 0.98)` - just below top

---

## 8. Title

```python
plt.suptitle(r'Resolution Study: $v^t_L = v^t_R = 0.9$, $t = 0.36$' + '\n',
            fontsize=13, fontweight='bold', y=1.01)
```

- Includes LaTeX math for velocity notation
- Shows simulation time
- Position: Above figure (`y=1.01`)

---

## 9. Grid

```python
ax.grid(True, alpha=0.2)
```

- Light grid for reference
- Alpha: 0.2 (subtle)

---

## 10. Output

```python
plt.tight_layout()
plt.savefig(save_path, dpi=150, bbox_inches='tight')
```

- DPI: 150 (good balance of quality and file size)
- Tight bounding box to remove whitespace

---

## 11. Time Consistency

**Critical**: All panels must use the same simulation time:

```python
SNAP_NUM = 12  # All resolutions use same snapshot
t_snapshot = SNAP_NUM * DT_OUTPUT  # = 0.36 for DT_OUTPUT=0.03

# Apply to all rows
resolutions = [
    (dir1, hllc_dir1, "800 + 800", SNAP_NUM),
    (dir2, hllc_dir2, "1,600 + 1,600", SNAP_NUM),
    # ...
]

# Analytical solution at same time
rho_ana, pres_ana, vx_ana, vt_ana = compute_analytical_solution(x, t_snapshot)
```

---

## 12. Complete Style Summary

```python
# Figure
fig, axes = plt.subplots(4, 3, figsize=(12, 14))

# Colors
EXACT_COLORS = {'P': 'darkgreen', 'n': 'darkred', 'vx': 'darkblue', 'vt': 'purple'}
HLLC_COLORS = {'P': 'lime', 'n': 'orange', 'vx': 'deepskyblue', 'vt': 'magenta'}
ANALYTICAL_COLOR = 'black'

# Markers
EXACT_STYLE = {'marker': 'x', 'ms': 3, 'mew': 0.8, 'alpha': 0.9, 'linestyle': 'None'}
HLLC_STYLE = {'marker': 'o', 'ms': 4, 'mfc': 'none', 'mew': 0.8, 'alpha': 0.7, 'linestyle': 'None'}
ANALYTICAL_STYLE = {'linestyle': '-', 'lw': 1.5, 'alpha': 0.8}

# Axes
XLIM = (-0.5, 0.5)
YLIMS = {0: (-0.05, 1.2), 1: (-0.02, 0.42), 2: (0.4, 1.02)}

# Fonts
TITLE_FONT = {'fontsize': 13, 'fontweight': 'bold'}
LABEL_FONT = {'fontsize': 11, 'fontweight': 'bold'}
COLUMN_TITLE_FONT = {'fontsize': 11}
LEGEND_FONT = {'fontsize': 11}

# Grid
GRID_ALPHA = 0.2

# Output
DPI = 150
```

---

## 13. Example Usage Prompt

When requesting a similar plot, you can say:

> "Create a Kitajima-style resolution comparison plot with:
> - 4 rows x 3 columns grid (figsize 12x14)
> - Rows labeled with resolution (e.g., '800 + 800')
> - Columns: Pressure/Density, Normal Velocity, Tangent Velocity
> - Dark colors for primary data (darkgreen, darkred, darkblue, purple)
> - Bright colors for comparison data (lime, orange, deepskyblue, magenta)
> - Black solid lines for analytical solution
> - x markers for primary data, hollow o markers for comparison
> - All panels at same time step with analytical overlay
> - Top-centered legend with 3 columns
> - Y-axis: P/1000 and n/5 scaled, vx in [0, 0.42], vt in [0.4, 1.0]"
