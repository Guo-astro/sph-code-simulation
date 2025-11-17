"""
Tests for the relativistic Riemann solver.
"""

import pytest
import numpy as np
from relativistic_riemann import (
    RelativisiticRiemannSolver,
    test_case_sr_sod,
    test_case_blast_wave,
    test_case_relativistic_shock,
)


class TestSolver:
    """Tests for RelativisiticRiemannSolver"""
    
    def test_initialization(self):
        """Test solver initialization"""
        solver = RelativisiticRiemannSolver(gamma=1.4)
        assert solver.gamma == 1.4
        
    def test_set_initial_states(self):
        """Test setting initial states"""
        solver = RelativisiticRiemannSolver(gamma=1.4)
        solver.set_initial_states(
            pl=1.0, rhol=1.0, vell=0.0,
            pr=0.1, rhor=0.125, velr=0.0
        )
        
        # Check that states are set
        assert solver.pl == 1.0
        assert solver.rhol == 1.0
        assert solver.vell == 0.0
        assert solver.pr == 0.1
        assert solver.rhor == 0.125
        assert solver.velr == 0.0
        
        # Check derived quantities are computed
        assert solver.ul is not None
        assert solver.hl is not None
        assert solver.csl is not None
        assert solver.wl is not None
        
    def test_sod_shock_tube(self):
        """Test SR Sod shock tube solution"""
        solver = RelativisiticRiemannSolver(gamma=1.4)
        solver.set_initial_states(
            pl=1.0, rhol=1.0, vell=0.0,
            pr=0.1, rhor=0.125, velr=0.0
        )
        
        x, p, rho, vel, u = solver.solve(t=0.4, n=400)
        
        # Check output shapes
        assert len(x) == 400
        assert len(p) == 400
        assert len(rho) == 400
        assert len(vel) == 400
        assert len(u) == 400
        
        # Check physical constraints
        assert np.all(rho > 0), "Density must be positive"
        assert np.all(p > 0), "Pressure must be positive"
        assert np.all(np.abs(vel) < 1.0), "Velocity must be less than speed of light"
        
        # Check that discontinuity is resolved
        assert rho.min() < 0.2, "Should have low density region"
        assert rho.max() > 0.9, "Should have high density region"
        
    def test_relativistic_shock(self):
        """Test relativistic shock with high velocity"""
        test = test_case_relativistic_shock()
        
        solver = RelativisiticRiemannSolver(gamma=test['gamma'])
        solver.set_initial_states(
            test['pl'], test['rhol'], test['vell'],
            test['pr'], test['rhor'], test['velr']
        )
        
        x, p, rho, vel, u = solver.solve(t=0.3, n=200)
        
        # Check that solution respects relativistic constraint
        assert np.all(np.abs(vel) < 1.0), "Velocity must be less than c"
        
        # Check physical positivity
        assert np.all(rho > 0)
        assert np.all(p > 0)
        
    def test_conservation_at_initial_time(self):
        """Test that solution matches initial conditions at t~0"""
        solver = RelativisiticRiemannSolver(gamma=1.4)
        solver.set_initial_states(
            pl=1.0, rhol=1.0, vell=0.0,
            pr=0.1, rhor=0.125, velr=0.0
        )
        
        # Solve at very small time
        x, p, rho, vel, u = solver.solve(t=0.001, n=400)
        
        # Left half should be close to initial left state
        left_indices = x < 0.48
        assert np.allclose(rho[left_indices], 1.0, rtol=0.01)
        assert np.allclose(p[left_indices], 1.0, rtol=0.01)
        
        # Right half should be close to initial right state  
        right_indices = x > 0.52
        assert np.allclose(rho[right_indices], 0.125, rtol=0.01)
        assert np.allclose(p[right_indices], 0.1, rtol=0.01)


class TestTestCases:
    """Tests for predefined test cases"""
    
    def test_sod_case(self):
        """Test SR Sod test case"""
        test = test_case_sr_sod()
        assert test['gamma'] == 1.4
        assert test['pl'] == 1.0
        assert test['pr'] == 0.1
        
    def test_blast_wave_case(self):
        """Test blast wave test case"""
        test = test_case_blast_wave()
        assert test['gamma'] == 5.0/3.0
        assert test['pl'] == 1000.0
        
    def test_relativistic_shock_case(self):
        """Test relativistic shock case"""
        test = test_case_relativistic_shock()
        assert test['gamma'] == 4.0/3.0
        assert test['vell'] == 0.9  # High velocity


class TestPhysicalConstraints:
    """Tests for physical constraints"""
    
    def test_speed_of_light_limit(self):
        """Test that velocities never exceed speed of light"""
        test_cases = [
            test_case_sr_sod(),
            test_case_blast_wave(),
            test_case_relativistic_shock(),
        ]
        
        for test in test_cases:
            solver = RelativisiticRiemannSolver(gamma=test['gamma'])
            solver.set_initial_states(
                test['pl'], test['rhol'], test['vell'],
                test['pr'], test['rhor'], test['velr']
            )
            
            x, p, rho, vel, u = solver.solve(t=test['tmax'], n=200)
            
            assert np.all(np.abs(vel) < 1.0), \
                f"Velocity exceeds c in {test['name']}"
    
    def test_positive_density_pressure(self):
        """Test that density and pressure remain positive"""
        solver = RelativisiticRiemannSolver(gamma=1.4)
        solver.set_initial_states(
            pl=1.0, rhol=1.0, vell=0.0,
            pr=0.1, rhor=0.125, velr=0.0
        )
        
        x, p, rho, vel, u = solver.solve(t=0.4, n=400)
        
        assert np.all(rho > 0), "Density must be positive"
        assert np.all(p > 0), "Pressure must be positive"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
