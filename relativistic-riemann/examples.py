#!/usr/bin/env python3
"""
Example usage of the relativistic Riemann solver package.
"""

from relativistic_riemann import (
    RelativisiticRiemannSolver,
    test_case_sr_sod,
    test_case_relativistic_shock,
)
import matplotlib.pyplot as plt


def example_basic_usage():
    """Basic example: solve and plot SR Sod problem"""
    print("Example 1: Basic SR Sod Shock Tube")
    print("=" * 50)
    
    # Create solver
    solver = RelativisiticRiemannSolver(gamma=1.4)
    solver.set_initial_states(
        pl=1.0, rhol=1.0, vell=0.0,
        pr=0.1, rhor=0.125, velr=0.0
    )
    
    # Solve at t=0.4
    x, p, rho, vel, u = solver.solve(t=0.4, n=400)
    
    # Plot results
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle('SR Sod Shock Tube at t=0.4', fontsize=14, weight='bold')
    
    axes[0, 0].plot(x, rho, 'b-', linewidth=2)
    axes[0, 0].set_ylabel('Density ρ')
    axes[0, 0].grid(True, alpha=0.3)
    
    axes[0, 1].plot(x, p, 'r-', linewidth=2)
    axes[0, 1].set_ylabel('Pressure p')
    axes[0, 1].grid(True, alpha=0.3)
    
    axes[1, 0].plot(x, vel, 'g-', linewidth=2)
    axes[1, 0].set_ylabel('Velocity v')
    axes[1, 0].set_xlabel('Position x')
    axes[1, 0].grid(True, alpha=0.3)
    
    axes[1, 1].plot(x, u, 'm-', linewidth=2)
    axes[1, 1].set_ylabel('Internal Energy u')
    axes[1, 1].set_xlabel('Position x')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('example_sod.png', dpi=150)
    print("Saved plot to example_sod.png\n")


def example_using_test_case():
    """Example using predefined test cases"""
    print("Example 2: Using Predefined Test Cases")
    print("=" * 50)
    
    # Load test case
    test = test_case_sr_sod()
    print(f"Test: {test['name']}")
    print(f"Gamma: {test['gamma']}")
    print(f"Left state: p={test['pl']}, ρ={test['rhol']}, v={test['vell']}")
    print(f"Right state: p={test['pr']}, ρ={test['rhor']}, v={test['velr']}")
    
    # Solve
    solver = RelativisiticRiemannSolver(gamma=test['gamma'])
    solver.set_initial_states(
        test['pl'], test['rhol'], test['vell'],
        test['pr'], test['rhor'], test['velr']
    )
    
    x, p, rho, vel, u = solver.solve(t=test['tmax'])
    
    print(f"Solution computed at t={test['tmax']}")
    print(f"Density range: [{rho.min():.4f}, {rho.max():.4f}]")
    print(f"Pressure range: [{p.min():.4f}, {p.max():.4f}]")
    print(f"Velocity range: [{vel.min():.4f}, {vel.max():.4f}]\n")


def example_relativistic_shock():
    """Example with high-velocity relativistic shock"""
    print("Example 3: Relativistic Shock (v=0.9c)")
    print("=" * 50)
    
    test = test_case_relativistic_shock()
    
    solver = RelativisiticRiemannSolver(gamma=test['gamma'])
    solver.set_initial_states(
        test['pl'], test['rhol'], test['vell'],
        test['pr'], test['rhor'], test['velr']
    )
    
    # Solve at multiple times
    times = [0.1, 0.2, 0.3, 0.5]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    colors = plt.cm.viridis([0.2, 0.4, 0.6, 0.9])
    
    for t, color in zip(times, colors):
        x, p, rho, vel, u = solver.solve(t=t)
        ax.plot(x, vel, label=f't={t:.2f}', color=color, linewidth=2)
    
    ax.set_xlabel('Position x', fontsize=12)
    ax.set_ylabel('Velocity v', fontsize=12)
    ax.set_title('Relativistic Shock Evolution (initial v=0.9c)', fontsize=14, weight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend()
    
    plt.tight_layout()
    plt.savefig('example_relativistic.png', dpi=150)
    print("Saved plot to example_relativistic.png\n")


def example_save_to_file():
    """Example saving solution to file"""
    print("Example 4: Saving Solution to File")
    print("=" * 50)
    
    solver = RelativisiticRiemannSolver(gamma=1.4)
    solver.set_initial_states(
        pl=1.0, rhol=1.0, vell=0.0,
        pr=0.1, rhor=0.125, velr=0.0
    )
    
    t = 0.4
    x, p, rho, vel, u = solver.solve(t=t)
    
    # Save to file
    solver.save_solution('solution.dat', x, p, rho, vel, u, t)
    print("Solution saved to solution.dat")
    print("File format: position, pressure, density, velocity, internal_energy\n")


if __name__ == "__main__":
    example_basic_usage()
    example_using_test_case()
    example_relativistic_shock()
    example_save_to_file()
    
    print("=" * 50)
    print("All examples completed successfully!")
