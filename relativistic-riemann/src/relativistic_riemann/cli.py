"""
Command-line interface for the relativistic Riemann solver.
"""

import sys
from .solver import RelativisiticRiemannSolver


def main():
    """Main CLI function - interactive mode"""
    print("Relativistic Riemann Problem Solver")
    print("=" * 50)
    
    # Get input parameters
    try:
        gamma = float(input("Adiabatic index of the gas: "))
        t = float(input("Time for the solution: "))
        
        print("\n-- LEFT STATE --")
        pl = float(input("  Pressure: "))
        rhol = float(input("  Density: "))
        vell = float(input("  Flow velocity: "))
        
        print("\n-- RIGHT STATE --")
        pr = float(input("  Pressure: "))
        rhor = float(input("  Density: "))
        velr = float(input("  Flow velocity: "))
    except (ValueError, EOFError, KeyboardInterrupt) as e:
        print(f"\nError reading input: {e}")
        sys.exit(1)
    
    # Create solver and solve
    solver = RelativisiticRiemannSolver(gamma)
    solver.set_initial_states(pl, rhol, vell, pr, rhor, velr)
    
    print("\nSolving...")
    try:
        x, pressure, density, velocity, internal_energy = solver.solve(t, n=400)
    except Exception as e:
        print(f"Error during solving: {e}")
        sys.exit(1)
    
    # Save solution
    output_file = 'solution.dat'
    solver.save_solution(output_file, x, pressure, density, velocity, internal_energy, t)
    
    print(f"Solution saved to {output_file}")
    

if __name__ == "__main__":
    main()
