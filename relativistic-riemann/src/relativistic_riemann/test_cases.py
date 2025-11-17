"""
Predefined test cases for relativistic Riemann problems.
"""


def test_case_sr_sod():
    """Standard SR Sod shock tube test"""
    return {
        'name': 'SR Sod Shock Tube',
        'gamma': 1.4,
        'pl': 1.0,
        'rhol': 1.0,
        'vell': 0.0,
        'pr': 0.1,
        'rhor': 0.125,
        'velr': 0.0,
        'tmax': 0.4
    }


def test_case_blast_wave():
    """Strong blast wave test"""
    return {
        'name': 'Blast Wave',
        'gamma': 5.0/3.0,
        'pl': 1000.0,
        'rhol': 1.0,
        'vell': 0.0,
        'pr': 0.01,
        'rhor': 1.0,
        'velr': 0.0,
        'tmax': 0.4
    }


def test_case_relativistic_shock():
    """Relativistic shock test with high velocities"""
    return {
        'name': 'Relativistic Shock',
        'gamma': 4.0/3.0,
        'pl': 10.0,
        'rhol': 1.0,
        'vell': 0.9,  # 0.9c
        'pr': 1.0,
        'rhor': 1.0,
        'velr': 0.0,
        'tmax': 0.5
    }


def test_case_two_shocks():
    """Two shock collision"""
    return {
        'name': 'Two Shocks Collision',
        'gamma': 1.4,
        'pl': 1.0,
        'rhol': 1.0,
        'vell': 0.5,
        'pr': 1.0,
        'rhor': 1.0,
        'velr': -0.5,
        'tmax': 0.6
    }
