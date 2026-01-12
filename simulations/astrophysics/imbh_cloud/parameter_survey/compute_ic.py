#!/usr/bin/env python3
"""
Compute initial conditions for parabolic IMBH-cloud encounter simulations.

For a parabolic orbit (E=0) with pericenter distance r_p:
- Angular momentum: L = sqrt(2 * G * M_BH * r_p)
- Orbital equation: r = p / (1 + cos(theta)), where p = 2 * r_p
- Velocity at distance r: v = sqrt(2 * G * M_BH / r)

Usage:
    python compute_ic.py --r_peri 0.4 --start_distance 15.0 --output config.json
"""

import argparse
import json
import numpy as np

# Physical constants in code units (L=pc, M=M_sun, V=km/s)
G = 0.00430091
M_BH = 100000.0  # 10^5 M_sun


def compute_parabolic_ic(r_peri: float, start_distance: float = 15.0) -> dict:
    """
    Compute initial conditions for a parabolic orbit.

    Args:
        r_peri: Pericenter distance in pc
        start_distance: Approximate starting distance from BH in pc

    Returns:
        Dictionary with position, velocity, and orbital parameters
    """
    # Semi-latus rectum for parabolic orbit
    p = 2 * r_peri

    # Angular momentum
    L = np.sqrt(2 * G * M_BH * r_peri)

    # Find true anomaly at starting distance
    # r = p / (1 + cos(theta))
    # For approaching orbit, use theta < 0 (before pericenter)

    # Solve for theta given r ~ start_distance
    # cos(theta) = p/r - 1
    cos_theta = p / start_distance - 1

    # Clamp to valid range
    if cos_theta < -1:
        cos_theta = -1
        print(f"Warning: start_distance too large, using theta = 180 deg")
    elif cos_theta > 1:
        cos_theta = 1
        print(f"Warning: start_distance too small, using theta = 0 deg")

    theta = np.arccos(cos_theta)

    # Actual distance at this theta
    r = p / (1 + np.cos(theta))

    # Use negative theta for approaching orbit (before pericenter)
    theta_start = -theta

    # Position in Cartesian (pericenter along +x axis)
    x = r * np.cos(theta_start)
    y = r * np.sin(theta_start)
    z = 0.0

    # Velocity components (orbital velocity for parabola)
    # v_r = sqrt(G*M/p) * sin(theta)
    # v_theta = sqrt(G*M/p) * (1 + cos(theta))
    v_factor = np.sqrt(G * M_BH / p)
    v_r = v_factor * np.sin(theta_start)
    v_theta = v_factor * (1 + np.cos(theta_start))

    # Convert to Cartesian
    # r_hat = (cos(theta), sin(theta), 0)
    # theta_hat = (-sin(theta), cos(theta), 0)
    vx = v_r * np.cos(theta_start) - v_theta * np.sin(theta_start)
    vy = v_r * np.sin(theta_start) + v_theta * np.cos(theta_start)
    vz = 0.0

    # Velocity at pericenter
    v_peri = np.sqrt(2 * G * M_BH / r_peri)

    # Total velocity (should equal sqrt(2GM/r))
    v_total = np.sqrt(vx**2 + vy**2)
    v_expected = np.sqrt(2 * G * M_BH / r)

    # Verify energy and angular momentum
    pos = np.array([x, y, z])
    vel = np.array([vx, vy, vz])
    L_calc = np.linalg.norm(np.cross(pos, vel))
    E_calc = 0.5 * np.dot(vel, vel) - G * M_BH / np.linalg.norm(pos)

    result = {
        "r_peri": r_peri,
        "p": p,
        "L": L,
        "theta_start_deg": np.degrees(theta_start),
        "position": [x, y, z],
        "velocity": [vx, vy, vz],
        "start_distance": r,
        "start_velocity": v_total,
        "v_pericenter": v_peri,
        "verification": {
            "L_calculated": L_calc,
            "L_expected": L,
            "E_calculated": E_calc,
            "E_expected": 0.0,
            "v_calculated": v_total,
            "v_expected": v_expected
        }
    }

    return result


def create_config_dict(ic: dict, name: str, cooling: bool = False) -> dict:
    """Create a full simulation config dictionary."""

    config = {
        "_comment": f"Parameter Survey: Parabolic orbit r_p={ic['r_peri']}pc - {'COOLING' if cooling else 'ADIABATIC'}",
        "_orbit": f"Parabolic (E=0), pericenter={ic['r_peri']}pc, v_peri={ic['v_pericenter']:.1f}km/s",
        "_units": "CODE units: L=pc, M=Msun, V=km/s",

        "_calculation": {
            "orbit_type": "parabolic (E=0)",
            "r_pericenter": f"{ic['r_peri']} pc",
            "v_pericenter": f"{ic['v_pericenter']:.2f} km/s",
            "start_distance": f"{ic['start_distance']:.2f} pc",
            "start_velocity": f"{ic['start_velocity']:.2f} km/s",
            "angular_momentum": f"{ic['L']:.4f}",
            "M_BH": "10^5 M_sun"
        },

        "name": name,
        "dimension": 3,

        "sample": "isothermal_bonnor_ebert",
        "resumeFromSnapshot": "simulations/astrophysics/imbh_cloud/results/phase2.5_compact_cooling/snapshot_0002.csv",
        "resetTimeOnResume": True,
        "useSnapshotSSOT": False,
        "N": 25,

        "M_cloud": 127.5,
        "T_cloud": 20.0,
        "n_center": 500.0,
        "xi_s": 3.0,
        "R_cloud": 1.55,
        "mu": 1.27,

        "gamma": 1.6667,
        "neighborNumber": 50,

        "G": G,
        "useGravity": True,

        "imbh_parameters": {
            "enabled": True,
            "use_code_units": True,
            "M_BH": M_BH,
            "BH_initial_position": {"0": 0.0, "1": 0.0, "2": 0.0},
            "BH_initial_velocity": {"0": 0.0, "1": 0.0, "2": 0.0},
            "softening_epsilon": 0.005,
            "sink_radius": 0.01,
            "enable_sink": True,
            "cloud_initial_position": {
                "0": ic["position"][0],
                "1": ic["position"][1],
                "2": ic["position"][2]
            },
            "cloud_initial_velocity": {
                "0": ic["velocity"][0],
                "1": ic["velocity"][1],
                "2": ic["velocity"][2]
            }
        },

        "sphAlgorithm": "gsph",
        "kernel": "wendland",
        "riemannSolver": "hll",
        "useGradH": False,

        "useBalsaraSwitch": True,
        "alpha": 1.0,

        "endTime": 6.0,
        "outputTime": 0.02,
        "energyOutputTime": 0.02,
        "cfl": 0.2,

        "useRelaxation": False,
        "useCooling": cooling,
        "outputFormat": "csv",
        "saveEnergyFile": True
    }

    if cooling:
        config["coolingFunction"] = "koyama_inutsuka_2000"

    return config


def main():
    parser = argparse.ArgumentParser(description="Compute IC for parabolic IMBH-cloud orbits")
    parser.add_argument("--r_peri", type=float, required=True, help="Pericenter distance in pc")
    parser.add_argument("--start_distance", type=float, default=15.0, help="Starting distance in pc")
    parser.add_argument("--output", type=str, help="Output JSON file path")
    parser.add_argument("--cooling", action="store_true", help="Generate cooling config")
    parser.add_argument("--name", type=str, help="Simulation name")

    args = parser.parse_args()

    # Compute initial conditions
    ic = compute_parabolic_ic(args.r_peri, args.start_distance)

    print("\n" + "="*60)
    print(f"Initial Conditions for r_p = {args.r_peri} pc")
    print("="*60)
    print(f"  Pericenter distance: {ic['r_peri']:.4f} pc")
    print(f"  Semi-latus rectum:   {ic['p']:.4f} pc")
    print(f"  Angular momentum:    {ic['L']:.4f}")
    print(f"  True anomaly:        {ic['theta_start_deg']:.2f} deg")
    print()
    print(f"  Starting position:   [{ic['position'][0]:.4f}, {ic['position'][1]:.4f}, {ic['position'][2]:.4f}] pc")
    print(f"  Starting velocity:   [{ic['velocity'][0]:.4f}, {ic['velocity'][1]:.4f}, {ic['velocity'][2]:.4f}] km/s")
    print(f"  |position|:          {ic['start_distance']:.4f} pc")
    print(f"  |velocity|:          {ic['start_velocity']:.4f} km/s")
    print()
    print(f"  Velocity at peri:    {ic['v_pericenter']:.2f} km/s")
    print()
    print("Verification:")
    print(f"  L calc/expected:     {ic['verification']['L_calculated']:.4f} / {ic['verification']['L_expected']:.4f}")
    print(f"  E (should be ~0):    {ic['verification']['E_calculated']:.8f}")
    print(f"  v calc/expected:     {ic['verification']['v_calculated']:.4f} / {ic['verification']['v_expected']:.4f}")

    if args.output:
        name = args.name or f"Parabolic r_p={args.r_peri}pc {'Cooling' if args.cooling else 'Adiabatic'}"
        config = create_config_dict(ic, name, args.cooling)

        with open(args.output, 'w') as f:
            json.dump(config, f, indent=4)
        print(f"\nConfig written to: {args.output}")


if __name__ == "__main__":
    main()
