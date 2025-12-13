#!/usr/bin/env python3
"""
Radice et al. (2018) Data Loader

Loads dynamical ejecta data from the Zenodo dataset:
"Binary Neutron Star Mergers: Mass Ejection, Electromagnetic Counterparts, 
and Nucleosynthesis"

Dataset DOI: 10.5281/zenodo.1298474
"""

import os
from pathlib import Path
from typing import Dict, Optional, Tuple, List
from dataclasses import dataclass

import numpy as np

try:
    import h5py
    HDF5_AVAILABLE = True
except ImportError:
    HDF5_AVAILABLE = False
    print("Warning: h5py not available, HDF5 loading disabled")


@dataclass
class EjectaProfile:
    """Container for ejecta profile data."""
    # Velocity distribution (1D)
    v_inf: np.ndarray          # Asymptotic velocity bins [c]
    dM_dv: np.ndarray          # Differential mass dM/dv [M_sun/c]
    
    # Angular distribution (1D)
    cos_theta: np.ndarray      # cos(theta) bins
    M_theta: np.ndarray        # Mass per angular bin [M_sun]
    v_theta: np.ndarray        # Average velocity per angle [c]
    
    # 2D velocity-angle distribution (optional)
    dM_dv_dOmega: Optional[np.ndarray] = None  # dM/(dv dOmega) [M_sun/c/sr]
    
    # Integrated quantities
    total_mass: float = 0.0    # Total ejecta mass [M_sun]
    avg_velocity: float = 0.0  # Mass-weighted average velocity [c]
    kinetic_energy: float = 0.0  # Total kinetic energy [erg]
    
    # Metadata
    model_name: str = ""
    eos: str = ""


class RadiceDataLoader:
    """
    Load and process Radice+2018 dynamical ejecta data.
    
    Usage:
        loader = RadiceDataLoader('/path/to/radice2018/data')
        profile = loader.load_model('SFHo_M140140')
    """
    
    # Physical constants
    M_SUN_CGS = 1.989e33  # g
    C_CGS = 2.998e10      # cm/s
    
    def __init__(self, data_dir: str):
        """
        Initialize loader with path to Radice+2018 data directory.
        
        Args:
            data_dir: Path to directory containing model subdirectories
        """
        self.data_dir = Path(data_dir)
        if not self.data_dir.exists():
            raise FileNotFoundError(f"Data directory not found: {data_dir}")
        
        self._available_models = self._scan_models()
    
    def _scan_models(self) -> List[str]:
        """Scan directory for available models."""
        models = []
        for item in self.data_dir.iterdir():
            if item.is_dir() and not item.name.startswith('.'):
                # Check for required files
                if (item / 'profile.txt').exists() or (item / 'hist_vinf.dat').exists():
                    models.append(item.name)
        return sorted(models)
    
    @property
    def available_models(self) -> List[str]:
        """List of available model names."""
        return self._available_models
    
    def load_model(self, model_name: str) -> EjectaProfile:
        """
        Load complete ejecta profile for a model.
        
        Args:
            model_name: Name of the model directory
            
        Returns:
            EjectaProfile with velocity and angular distributions
        """
        model_dir = self.data_dir / model_name
        if not model_dir.exists():
            raise ValueError(f"Model not found: {model_name}. "
                           f"Available: {self.available_models}")
        
        # Load 1D velocity distribution
        v_inf, dM_dv = self._load_hist_vinf(model_dir)
        
        # Load angular profile
        cos_theta, M_theta, v_theta = self._load_profile(model_dir)
        
        # Load 2D distribution if available
        dM_dv_dOmega = None
        if HDF5_AVAILABLE:
            dM_dv_dOmega = self._load_hist_vinf_theta(model_dir, v_inf, cos_theta)
        
        # Compute integrated quantities
        total_mass = np.trapz(dM_dv, v_inf)
        avg_velocity = np.trapz(v_inf * dM_dv, v_inf) / total_mass if total_mass > 0 else 0
        
        # Kinetic energy: E_kin = (Γ-1) M c²
        gamma = 1.0 / np.sqrt(1.0 - v_inf**2)
        dE_dv = (gamma - 1.0) * dM_dv * self.M_SUN_CGS * self.C_CGS**2
        kinetic_energy = np.trapz(dE_dv, v_inf)
        
        # Parse EOS from model name
        eos = model_name.split('_')[0] if '_' in model_name else 'unknown'
        
        return EjectaProfile(
            v_inf=v_inf,
            dM_dv=dM_dv,
            cos_theta=cos_theta,
            M_theta=M_theta,
            v_theta=v_theta,
            dM_dv_dOmega=dM_dv_dOmega,
            total_mass=total_mass,
            avg_velocity=avg_velocity,
            kinetic_energy=kinetic_energy,
            model_name=model_name,
            eos=eos,
        )
    
    def _load_hist_vinf(self, model_dir: Path) -> Tuple[np.ndarray, np.ndarray]:
        """Load velocity histogram from hist_vinf.dat."""
        filepath = model_dir / 'hist_vinf.dat'
        if not filepath.exists():
            # Return empty arrays if file doesn't exist
            return np.array([0.1, 0.2, 0.3]), np.array([0.0, 0.0, 0.0])
        
        data = np.loadtxt(filepath)
        v_inf = data[:, 0]   # velocity [c]
        dM_dv = data[:, 1]   # dM/dv [M_sun/c]
        
        return v_inf, dM_dv
    
    def _load_profile(self, model_dir: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Load angular profile from profile.txt."""
        filepath = model_dir / 'profile.txt'
        if not filepath.exists():
            # Return uniform angular distribution
            cos_theta = np.linspace(-1, 1, 20)
            M_theta = np.ones_like(cos_theta) * 0.001
            v_theta = np.ones_like(cos_theta) * 0.2
            return cos_theta, M_theta, v_theta
        
        data = np.loadtxt(filepath)
        cos_theta = data[:, 0]  # cos(theta)
        M_theta = data[:, 1]    # mass [M_sun]
        v_theta = data[:, 2]    # velocity [c]
        
        return cos_theta, M_theta, v_theta
    
    def _load_hist_vinf_theta(self, model_dir: Path, 
                              v_inf: np.ndarray, 
                              cos_theta: np.ndarray) -> Optional[np.ndarray]:
        """Load 2D velocity-angle histogram from HDF5."""
        filepath = model_dir / 'hist_vinf_theta.h5'
        if not filepath.exists() or not HDF5_AVAILABLE:
            return None
        
        try:
            with h5py.File(filepath, 'r') as f:
                # Check available keys
                if 'dM_dv_dOmega' in f:
                    return f['dM_dv_dOmega'][:]
                elif 'mass' in f:
                    return f['mass'][:]
                else:
                    print(f"Warning: Unknown HDF5 structure in {filepath}")
                    return None
        except Exception as e:
            print(f"Warning: Could not load HDF5 file: {e}")
            return None
    
    def get_spherical_average(self, profile: EjectaProfile) -> Tuple[np.ndarray, np.ndarray]:
        """
        Compute spherically averaged mass-velocity distribution.
        
        Returns:
            v: Velocity bins [c]
            M_v: Mass per velocity bin [M_sun]
        """
        # Simply return the 1D velocity histogram
        return profile.v_inf, profile.dM_dv * np.gradient(profile.v_inf)
    
    def print_summary(self, profile: EjectaProfile):
        """Print summary of ejecta profile."""
        print(f"=== Ejecta Profile: {profile.model_name} ===")
        print(f"EOS: {profile.eos}")
        print(f"Total mass: {profile.total_mass:.4f} M_sun")
        print(f"Average velocity: {profile.avg_velocity:.3f} c")
        print(f"Kinetic energy: {profile.kinetic_energy:.3e} erg")
        print(f"Velocity range: {profile.v_inf.min():.3f} - {profile.v_inf.max():.3f} c")
        print(f"Angular bins: {len(profile.cos_theta)}")
        if profile.dM_dv_dOmega is not None:
            print(f"2D distribution: {profile.dM_dv_dOmega.shape}")


def main():
    """Test the data loader."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Load Radice+2018 ejecta data')
    parser.add_argument('data_dir', help='Path to Radice+2018 data directory')
    parser.add_argument('--model', '-m', help='Specific model to load')
    parser.add_argument('--list', '-l', action='store_true', help='List available models')
    
    args = parser.parse_args()
    
    loader = RadiceDataLoader(args.data_dir)
    
    if args.list:
        print("Available models:")
        for model in loader.available_models:
            print(f"  {model}")
        return
    
    if args.model:
        profile = loader.load_model(args.model)
        loader.print_summary(profile)
    else:
        print(f"Found {len(loader.available_models)} models")
        print("Use --model NAME to load a specific model")
        print("Use --list to see available models")


if __name__ == '__main__':
    main()
