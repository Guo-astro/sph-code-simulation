#!/usr/bin/env python3
"""
Generate comparison plots: Exact vs HLLC solver with analytical solution overlay
Uses the SR Riemann solver from src/srgsph/reference/riemann_solver.py
"""

import sys
sys.path.insert(0, '/Users/guo/Downloads/sphcode/src/srgsph/reference')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import glob
from scipy.optimize import brentq
from scipy.integrate import solve_ivp

# Configuration
RESULTS_DIR = Path("/Users/guo/Downloads/sphcode/simulations/relativistic/sr_tangent/results/kitajima_paper_resolution")
DT_OUTPUT = 0.03

# Initial conditions for tangent velocity test
P_L, RHO_L, VX_L, VT_L = 1000.0, 1.0, 0.0, 0.9
P_R, RHO_R, VX_R, VT_R = 0.01, 1.0, 0.0, 0.9
GAMMA = 5.0/3.0

class SRRiemannSolver:
    """Special Relativistic Riemann Problem Solver with spatial sampling"""

    def __init__(self, rho_L, P_L, vx_L, vt_L, rho_R, P_R, vx_R, vt_R, gamma=5.0/3.0):
        self.gamma = gamma

        # Left state
        self.rho_L = rho_L
        self.P_L = P_L
        self.vx_L = vx_L
        self.vt_L = vt_L
        self.h_L = self._enthalpy(P_L, rho_L)
        self.cs_L = self._sound_speed(P_L, rho_L)
        self.W_L = 1.0 / np.sqrt(1.0 - vx_L**2 - vt_L**2)

        # Right state
        self.rho_R = rho_R
        self.P_R = P_R
        self.vx_R = vx_R
        self.vt_R = vt_R
        self.h_R = self._enthalpy(P_R, rho_R)
        self.cs_R = self._sound_speed(P_R, rho_R)
        self.W_R = 1.0 / np.sqrt(1.0 - vx_R**2 - vt_R**2)

        # Solve for star state
        self._solve_star_state()

    def _enthalpy(self, P, rho):
        """Specific enthalpy h = 1 + gamma/(gamma-1) * P/rho"""
        return 1.0 + (self.gamma / (self.gamma - 1.0)) * (P / rho)

    def _sound_speed(self, P, rho):
        """Sound speed cs = sqrt(gamma * P / (rho * h))"""
        h = self._enthalpy(P, rho)
        return np.sqrt(self.gamma * P / (rho * h))

    def _solve_star_state(self):
        """Find P* and v* by matching wave curves"""

        def velocity_jump_left(P_star):
            """Velocity behind left wave (shock or rarefaction)"""
            if P_star > self.P_L:
                return self._shock_velocity(P_star, self.P_L, self.rho_L, self.h_L,
                                           self.vx_L, self.vt_L, self.W_L, 'left')
            else:
                return self._rarefaction_velocity(P_star, self.P_L, self.rho_L,
                                                 self.vx_L, self.vt_L, 'left')

        def velocity_jump_right(P_star):
            """Velocity behind right wave (shock or rarefaction)"""
            if P_star > self.P_R:
                return self._shock_velocity(P_star, self.P_R, self.rho_R, self.h_R,
                                           self.vx_R, self.vt_R, self.W_R, 'right')
            else:
                return self._rarefaction_velocity(P_star, self.P_R, self.rho_R,
                                                 self.vx_R, self.vt_R, 'right')

        def func(P):
            vL = velocity_jump_left(P)
            vR = velocity_jump_right(P)
            return vL - vR

        # Bracket the root
        P_min = min(self.P_L, self.P_R) * 1e-6
        P_max = max(self.P_L, self.P_R) * 10.0

        try:
            self.P_star = brentq(func, P_min, P_max, xtol=1e-10)
        except:
            self.P_star = 0.5 * (self.P_L + self.P_R)

        self.vx_star = velocity_jump_left(self.P_star)

        # Compute star state densities
        if self.P_star > self.P_L:
            self.rho_star_L = self._shock_density(self.P_star, self.P_L, self.rho_L, self.h_L)
        else:
            self.rho_star_L = self.rho_L * (self.P_star / self.P_L)**(1.0/self.gamma)

        if self.P_star > self.P_R:
            self.rho_star_R = self._shock_density(self.P_star, self.P_R, self.rho_R, self.h_R)
        else:
            self.rho_star_R = self.rho_R * (self.P_star / self.P_R)**(1.0/self.gamma)

        # Compute tangential velocities in star region
        self.h_star_L = self._enthalpy(self.P_star, self.rho_star_L)
        self.h_star_R = self._enthalpy(self.P_star, self.rho_star_R)

        # K-invariant: h * W * vt = const across waves
        K_L = self.h_L * self.W_L * self.vt_L
        K_R = self.h_R * self.W_R * self.vt_R

        W_star = 1.0 / np.sqrt(1.0 - self.vx_star**2)  # Approximate
        self.vt_star_L = K_L / (self.h_star_L * W_star) if self.h_star_L * W_star > 0 else self.vt_L
        self.vt_star_R = K_R / (self.h_star_R * W_star) if self.h_star_R * W_star > 0 else self.vt_R

        # Wave speeds
        self._compute_wave_speeds()

    def _shock_velocity(self, P_star, P_a, rho_a, h_a, vx_a, vt_a, W_a, side):
        """Post-shock velocity from Taub adiabat"""
        A = (self.gamma - 1.0) * (P_a - P_star) / (self.gamma * P_star)
        B = h_a * (P_a - P_star) / rho_a

        qa = 1.0 + A
        qb = -A
        qc = B - h_a**2

        delta = qb**2 - 4*qa*qc
        if delta < 0:
            return vx_a

        h_star = (-qb + np.sqrt(delta)) / (2*qa)
        rho_star = (self.gamma / (self.gamma - 1.0)) * P_star / (h_star - 1.0)

        denom = h_star/rho_star - h_a/rho_a
        if abs(denom) < 1e-15:
            return vx_a

        j2 = -(P_star - P_a) / denom
        if j2 <= 0:
            return vx_a

        j = np.sqrt(j2)
        D_a = rho_a * W_a

        term = j**2 + D_a**2 * (1.0 - vx_a**2)
        sqrt_term = np.sqrt(term)
        denom_vs = D_a**2 + j**2

        if side == 'left':
            Vs = (D_a**2 * vx_a - j * sqrt_term) / denom_vs
        else:
            Vs = (D_a**2 * vx_a + j * sqrt_term) / denom_vs

        Ws = 1.0 / np.sqrt(1.0 - Vs**2)
        j_signed = j if side == 'right' else -j

        num = h_a * W_a * vx_a + Ws * (P_star - P_a) / j_signed
        den = h_a * W_a + (P_star - P_a) * (Ws * vx_a / j_signed + 1.0 / D_a)

        if abs(den) < 1e-15:
            return vx_a

        return num / den

    def _shock_density(self, P_star, P_a, rho_a, h_a):
        """Post-shock density"""
        A = (self.gamma - 1.0) * (P_a - P_star) / (self.gamma * P_star)
        B = h_a * (P_a - P_star) / rho_a

        qa = 1.0 + A
        qb = -A
        qc = B - h_a**2

        delta = qb**2 - 4*qa*qc
        if delta < 0:
            return rho_a

        h_star = (-qb + np.sqrt(delta)) / (2*qa)
        return (self.gamma / (self.gamma - 1.0)) * P_star / (h_star - 1.0)

    def _rarefaction_velocity(self, P_star, P_a, rho_a, vx_a, vt_a, side):
        """Velocity change across rarefaction using ODE integration"""
        if abs(vt_a) < 1e-10:
            # Analytical for vt = 0
            return self._rarefaction_velocity_analytical(P_star, P_a, rho_a, vx_a, side)

        # K-invariant
        h_a = self._enthalpy(P_a, rho_a)
        W_a = 1.0 / np.sqrt(1.0 - vx_a**2 - vt_a**2)
        Kt = h_a * W_a * vt_a

        const_entropy = P_a / (rho_a**self.gamma)
        sign = -1.0 if side == 'left' else 1.0

        def deriv(P, y):
            vx = y[0]
            rho = (P / const_entropy)**(1.0/self.gamma)
            h = self._enthalpy(P, rho)

            W2 = (1.0 + (Kt/h)**2) / (1.0 - vx**2)
            W = np.sqrt(W2)
            vt = Kt / (h * W)

            cs2 = (self.gamma * P) / (rho * h)
            cs = np.sqrt(cs2)

            v2 = vx**2 + vt**2
            term1 = vx * (1.0 - cs2)
            term2_inner = (1.0 - v2) * (1.0 - v2*cs2 - vx**2 * (1.0 - cs2))
            if term2_inner < 0: term2_inner = 0
            term2 = cs * np.sqrt(term2_inner)
            denom = 1.0 - v2 * cs2

            xi_sign = -1.0 if side == 'left' else 1.0
            xi = (term1 + xi_sign * term2) / denom

            g = (vt**2 * (xi**2 - 1.0)) / ((1.0 - xi * vx)**2 + 1e-20)

            dvdP = sign * (1.0 / (rho * h * W2 * cs + 1e-20)) * (1.0 / np.sqrt(1.0 + g))
            return [dvdP]

        try:
            res = solve_ivp(deriv, [P_a, P_star], [vx_a], rtol=1e-8, atol=1e-10)
            return res.y[0][-1]
        except:
            return vx_a

    def _rarefaction_velocity_analytical(self, P_star, P_a, rho_a, vx_a, side):
        """Analytical solution for rarefaction with vt=0"""
        const_entropy = P_a / (rho_a**self.gamma)
        rho_star = (P_star / const_entropy)**(1.0/self.gamma)

        h_a = self._enthalpy(P_a, rho_a)
        h_star = self._enthalpy(P_star, rho_star)

        cs_a = np.sqrt(self.gamma * P_a / (rho_a * h_a))
        cs_star = np.sqrt(self.gamma * P_star / (rho_star * h_star))

        sign = 1.0 if side == 'left' else -1.0
        sqrt_g = np.sqrt(self.gamma - 1.0)

        term_v = (1.0 + vx_a) / (1.0 - vx_a + 1e-20)
        term_ca = (sqrt_g + cs_a) / (sqrt_g - cs_a + 1e-10)
        term_cb = (sqrt_g + cs_star) / (sqrt_g - cs_star + 1e-10)

        exponent = sign * 2.0 / sqrt_g
        A = term_v * (term_ca / term_cb)**exponent

        return (A - 1.0) / (A + 1.0)

    def _compute_wave_speeds(self):
        """Compute characteristic wave speeds"""
        # Left wave
        if self.P_star > self.P_L:
            # Shock - use shock speed
            j2 = -(self.P_star - self.P_L) / (self.h_star_L/self.rho_star_L - self.h_L/self.rho_L + 1e-20)
            if j2 > 0:
                j = np.sqrt(j2)
                D_L = self.rho_L * self.W_L
                term = j**2 + D_L**2 * (1.0 - self.vx_L**2)
                self.S_L = (D_L**2 * self.vx_L - j * np.sqrt(term)) / (D_L**2 + j**2)
            else:
                self.S_L = self.vx_L - self.cs_L
        else:
            # Rarefaction
            self.xi_L_head = (self.vx_L - self.cs_L) / (1.0 - self.vx_L * self.cs_L)
            cs_star_L = self._sound_speed(self.P_star, self.rho_star_L)
            self.xi_L_tail = (self.vx_star - cs_star_L) / (1.0 - self.vx_star * cs_star_L)
            self.S_L = self.xi_L_head

        # Right wave
        if self.P_star > self.P_R:
            j2 = -(self.P_star - self.P_R) / (self.h_star_R/self.rho_star_R - self.h_R/self.rho_R + 1e-20)
            if j2 > 0:
                j = np.sqrt(j2)
                D_R = self.rho_R * self.W_R
                term = j**2 + D_R**2 * (1.0 - self.vx_R**2)
                self.S_R = (D_R**2 * self.vx_R + j * np.sqrt(term)) / (D_R**2 + j**2)
            else:
                self.S_R = self.vx_R + self.cs_R
        else:
            self.xi_R_head = (self.vx_R + self.cs_R) / (1.0 + self.vx_R * self.cs_R)
            cs_star_R = self._sound_speed(self.P_star, self.rho_star_R)
            self.xi_R_tail = (self.vx_star + cs_star_R) / (1.0 + self.vx_star * cs_star_R)
            self.S_R = self.xi_R_head

        # Contact speed
        self.S_contact = self.vx_star

    def sample(self, x, t):
        """Sample solution at position x and time t"""
        if t <= 0:
            if x < 0:
                return {'P': self.P_L, 'rho': self.rho_L, 'vx': self.vx_L, 'vt': self.vt_L}
            else:
                return {'P': self.P_R, 'rho': self.rho_R, 'vx': self.vx_R, 'vt': self.vt_R}

        xi = x / t  # Similarity variable

        # Left state
        if self.P_star > self.P_L:
            # Left shock
            if xi < self.S_L:
                return {'P': self.P_L, 'rho': self.rho_L, 'vx': self.vx_L, 'vt': self.vt_L}
        else:
            # Left rarefaction
            if xi < self.xi_L_head:
                return {'P': self.P_L, 'rho': self.rho_L, 'vx': self.vx_L, 'vt': self.vt_L}
            elif xi < self.xi_L_tail:
                return self._sample_rarefaction(xi, 'left')

        # Star region
        if xi < self.S_contact:
            return {'P': self.P_star, 'rho': self.rho_star_L, 'vx': self.vx_star, 'vt': self.vt_star_L}

        # Right side of contact
        if self.P_star > self.P_R:
            # Right shock
            if xi < self.S_R:
                return {'P': self.P_star, 'rho': self.rho_star_R, 'vx': self.vx_star, 'vt': self.vt_star_R}
        else:
            # Right rarefaction
            if xi < self.xi_R_tail:
                return {'P': self.P_star, 'rho': self.rho_star_R, 'vx': self.vx_star, 'vt': self.vt_star_R}
            elif xi < self.xi_R_head:
                return self._sample_rarefaction(xi, 'right')

        # Right state
        return {'P': self.P_R, 'rho': self.rho_R, 'vx': self.vx_R, 'vt': self.vt_R}

    def _sample_rarefaction(self, xi, side):
        """Sample inside rarefaction fan"""
        if side == 'left':
            P_a, rho_a, vx_a, vt_a = self.P_L, self.rho_L, self.vx_L, self.vt_L
            P_b = self.P_star
        else:
            P_a, rho_a, vx_a, vt_a = self.P_R, self.rho_R, self.vx_R, self.vt_R
            P_b = self.P_star

        # Binary search for P at this xi
        def xi_at_P(P):
            rho = rho_a * (P / P_a)**(1.0/self.gamma)
            vx = self._rarefaction_velocity(P, P_a, rho_a, vx_a, vt_a, side)
            cs = self._sound_speed(P, rho)
            if side == 'left':
                return (vx - cs) / (1.0 - vx * cs)
            else:
                return (vx + cs) / (1.0 + vx * cs)

        try:
            P = brentq(lambda P: xi_at_P(P) - xi, min(P_a, P_b), max(P_a, P_b))
        except:
            P = 0.5 * (P_a + P_b)

        rho = rho_a * (P / P_a)**(1.0/self.gamma)
        vx = self._rarefaction_velocity(P, P_a, rho_a, vx_a, vt_a, side)

        # K-invariant for vt
        h_a = self._enthalpy(P_a, rho_a)
        W_a = 1.0 / np.sqrt(1.0 - vx_a**2 - vt_a**2 + 1e-20)
        Kt = h_a * W_a * vt_a
        h = self._enthalpy(P, rho)
        W = 1.0 / np.sqrt(1.0 - vx**2 + 1e-20)
        vt = Kt / (h * W) if h * W > 0 else vt_a

        return {'P': P, 'rho': rho, 'vx': vx, 'vt': vt}


def get_analytical_solution(t, x_array):
    """Compute analytical solution"""
    try:
        rp = SRRiemannSolver(
            rho_L=RHO_L, P_L=P_L, vx_L=VX_L, vt_L=VT_L,
            rho_R=RHO_R, P_R=P_R, vx_R=VX_R, vt_R=VT_R,
            gamma=GAMMA
        )

        pres = np.zeros_like(x_array)
        dens = np.zeros_like(x_array)
        vx = np.zeros_like(x_array)
        vt = np.zeros_like(x_array)

        for i, x in enumerate(x_array):
            state = rp.sample(x, t)
            pres[i] = state['P']
            dens[i] = state['rho']
            vx[i] = state['vx']
            vt[i] = state['vt']

        return {'x': x_array, 'pres': pres, 'dens': dens, 'vx': vx, 'vt': vt}
    except Exception as e:
        print(f"Error computing analytical solution: {e}")
        return None


def load_snapshot(output_dir):
    """Load last snapshot from output directory."""
    files = sorted(glob.glob(str(output_dir / 'snapshot_*.csv')))
    if not files:
        return None, None
    f = files[-1]
    snap_num = int(Path(f).stem.split('_')[-1])
    t = snap_num * DT_OUTPUT
    df = pd.read_csv(f, comment='#').dropna(subset=['pos_x'])
    return df, t


def plot_solver_comparison(exact_dir, hllc_dir, N_left, N_right, save_path):
    """Create comparison plot with analytical solution."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    exact_df, t_exact = load_snapshot(exact_dir)
    hllc_df, t_hllc = load_snapshot(hllc_dir)

    t = t_exact if t_exact else t_hllc if t_hllc else 0.45

    # Analytical solution
    x_anal = np.linspace(-0.4, 0.4, 1000)
    anal = get_analytical_solution(t, x_anal)

    # Plot settings
    ms = 4
    alpha = 0.7

    # === Pressure ===
    ax = axes[0, 0]
    if anal is not None:
        ax.plot(anal['x'], anal['pres'] / 1000, 'k-', lw=2.5, label='Analytical', zorder=5)
    if exact_df is not None:
        ax.scatter(exact_df['pos_x'], exact_df['pres'] / 1000, s=ms, c='blue', alpha=alpha, label='Exact Riemann', zorder=2)
    if hllc_df is not None:
        ax.scatter(hllc_df['pos_x'], hllc_df['pres'] / 1000, s=ms, c='red', alpha=alpha, label='HLLC', zorder=3)
    ax.set_ylabel('P / 1000', fontsize=12)
    ax.set_xlim(-0.4, 0.4)
    ax.set_ylim(-0.05, 1.1)
    ax.set_title('Pressure', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10, loc='upper right', markerscale=3)
    ax.grid(True, alpha=0.3)

    # === Number Density ===
    ax = axes[0, 1]
    if anal is not None:
        ax.plot(anal['x'], anal['dens'] / 5, 'k-', lw=2.5, label='Analytical', zorder=5)
    if exact_df is not None:
        ax.scatter(exact_df['pos_x'], exact_df['dens'] / 5, s=ms, c='blue', alpha=alpha, label='Exact Riemann', zorder=2)
    if hllc_df is not None:
        ax.scatter(hllc_df['pos_x'], hllc_df['dens'] / 5, s=ms, c='red', alpha=alpha, label='HLLC', zorder=3)
    ax.set_ylabel('n / 5', fontsize=12)
    ax.set_xlim(-0.4, 0.4)
    ax.set_ylim(-0.05, 1.1)
    ax.set_title('Number Density', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10, loc='upper right', markerscale=3)
    ax.grid(True, alpha=0.3)

    # === Normal Velocity ===
    ax = axes[1, 0]
    if anal is not None:
        ax.plot(anal['x'], anal['vx'], 'k-', lw=2.5, label='Analytical', zorder=5)
    if exact_df is not None:
        ax.scatter(exact_df['pos_x'], exact_df['vel_x'], s=ms, c='blue', alpha=alpha, label='Exact Riemann', zorder=2)
    if hllc_df is not None:
        ax.scatter(hllc_df['pos_x'], hllc_df['vel_x'], s=ms, c='red', alpha=alpha, label='HLLC', zorder=3)
    ax.set_ylabel(r'$v^x$', fontsize=12)
    ax.set_xlabel('x', fontsize=12)
    ax.set_xlim(-0.4, 0.4)
    ax.set_ylim(-0.02, 0.45)
    ax.set_title('Normal Velocity', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10, loc='upper left', markerscale=3)
    ax.grid(True, alpha=0.3)

    # === Tangent Velocity ===
    ax = axes[1, 1]
    if anal is not None:
        ax.plot(anal['x'], anal['vt'], 'k-', lw=2.5, label='Analytical', zorder=5)
    if exact_df is not None and 'vel_t' in exact_df.columns:
        ax.scatter(exact_df['pos_x'], exact_df['vel_t'], s=ms, c='blue', alpha=alpha, label='Exact Riemann', zorder=2)
    if hllc_df is not None and 'vel_t' in hllc_df.columns:
        ax.scatter(hllc_df['pos_x'], hllc_df['vel_t'], s=ms, c='red', alpha=alpha, label='HLLC', zorder=3)
    ax.set_ylabel(r'$v^t$', fontsize=12)
    ax.set_xlabel('x', fontsize=12)
    ax.set_xlim(-0.4, 0.4)
    ax.set_ylim(0.4, 1.02)
    ax.set_title('Tangent Velocity', fontsize=12, fontweight='bold')
    ax.legend(fontsize=10, loc='lower left', markerscale=3)
    ax.grid(True, alpha=0.3)

    plt.suptitle(f'SR-GSPH: Exact vs HLLC Riemann Solver\n'
                 f'Tangent Velocity Test ($v^t_L = v^t_R = 0.9$), '
                 f'N_left={N_left}, N_right={N_right}, t={t:.3f}',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()


def plot_error_comparison(save_path):
    """Plot error vs resolution for both solvers."""
    V_SHOCK = 0.4451

    def find_shock_position(x, dens):
        sort_idx = np.argsort(x)
        x_sorted = x[sort_idx]
        dens_sorted = dens[sort_idx]
        right_mask = x_sorted > 0.02
        x_right = x_sorted[right_mask]
        dens_right = dens_sorted[right_mask]
        if len(dens_right) < 10:
            return None
        high_dens_mask = dens_right > 2.5
        if any(high_dens_mask):
            return x_right[high_dens_mask][-1]
        return None

    def get_error(output_dir):
        df, t = load_snapshot(output_dir)
        if df is None:
            return None
        x_shock = find_shock_position(df['pos_x'].values, df['dens'].values)
        if x_shock is None:
            return None
        V_measured = x_shock / t
        return (V_measured - V_SHOCK) / V_SHOCK * 100

    n_lefts = [1600, 3200, 6400, 12800]
    exact_errs = []
    hllc_errs = []

    for n_left in n_lefts:
        exact_errs.append(get_error(RESULTS_DIR / f'exact_{n_left}_1600'))
        hllc_errs.append(get_error(RESULTS_DIR / f'hllc_{n_left}_1600'))

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.semilogx(n_lefts, exact_errs, 'bo-', ms=12, lw=2.5, label='Exact Riemann Solver')
    ax.semilogx(n_lefts, hllc_errs, 'rs-', ms=12, lw=2.5, label='HLLC Solver')
    ax.axhline(y=0, color='green', ls='--', lw=2, label='Analytical (0% error)')

    for n, e_err, h_err in zip(n_lefts, exact_errs, hllc_errs):
        if e_err is not None:
            ax.annotate(f'{e_err:+.1f}%', (n, e_err), textcoords="offset points",
                       xytext=(-15, 10), ha='center', fontsize=10, color='blue')
        if h_err is not None:
            ax.annotate(f'{h_err:+.1f}%', (n, h_err), textcoords="offset points",
                       xytext=(15, -15), ha='center', fontsize=10, color='red')

    ax.set_xlabel('N_left (Rarefaction Side Particles)', fontsize=12)
    ax.set_ylabel('Shock Speed Error (%)', fontsize=12)
    ax.set_title('Kitajima Hypothesis: Left Resolution Effect on Shock Speed\n'
                'N_right fixed at 1600, $v^t_L = v^t_R = 0.9$', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11, loc='upper right')
    ax.grid(True, alpha=0.3)
    ax.set_xticks(n_lefts)
    ax.set_xticklabels([str(n) for n in n_lefts])

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Saved: {save_path}")
    plt.close()


if __name__ == '__main__':
    print("Generating comparison plots with analytical solution...")

    # Solver comparison at baseline
    plot_solver_comparison(
        RESULTS_DIR / "exact_1600_1600",
        RESULTS_DIR / "hllc_1600_1600",
        1600, 1600,
        RESULTS_DIR / 'solver_comparison_1600.png'
    )

    # Solver comparison at high resolution
    plot_solver_comparison(
        RESULTS_DIR / "exact_12800_1600",
        RESULTS_DIR / "hllc_12800_1600",
        12800, 1600,
        RESULTS_DIR / 'solver_comparison_12800.png'
    )

    # Error vs resolution
    plot_error_comparison(RESULTS_DIR / 'error_vs_resolution.png')

    print("\nAll plots generated!")
