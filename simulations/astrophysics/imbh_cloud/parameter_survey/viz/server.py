#!/usr/bin/env python3
"""
Auto-discovery HTTP server for IMBH-Cloud visualization.

This server automatically scans the parameter_survey directory for simulation
datasets on startup and whenever datasets.json is requested.

Usage:
    cd parameter_survey
    python viz/server.py [port]

    # Default port is 8888
    # Then open http://localhost:8888/viz/
"""
import http.server
import socketserver
import os
import sys
import json
import re
from pathlib import Path
from functools import partial

# Import scan logic
def extract_orbit_info(dirname):
    """Extract orbit info from directory name."""
    match = re.match(r'parabolic_rp(\d+\.?\d*)pc(?:_(.+))?', dirname)
    if match:
        r_peri = float(match.group(1))
        suffix = match.group(2) or ''
        return {'r_peri': r_peri, 'orbit_id': f'rp{r_peri}pc', 'suffix': suffix}
    return None

def read_config(config_path):
    """Read simulation config JSON file."""
    try:
        with open(config_path, 'r') as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return None

def count_snapshots(results_dir):
    """Count HDF5 snapshots in results directory."""
    if not results_dir.exists():
        return 0
    return len(list(results_dir.glob('snapshot_*.h5')))

def scan_datasets(survey_dir):
    """Scan parameter_survey directory and return datasets config."""
    skip_dirs = {'viz', '__pycache__', '.git'}
    all_datasets = []

    for item in sorted(survey_dir.iterdir()):
        if not item.is_dir() or item.name in skip_dirs or item.name.startswith('.'):
            continue

        orbit_info = extract_orbit_info(item.name)
        if not orbit_info:
            continue

        config_dir = item / 'config'

        for physics in ['adiabatic', 'cooling']:
            config_path = config_dir / f'{physics}.json'
            results_dir = item / physics / 'results'

            if not config_path.exists() and not results_dir.exists():
                continue

            config = read_config(config_path) or {}
            n_snapshots = count_snapshots(results_dir)
            if n_snapshots == 0:
                continue

            suffix = orbit_info['suffix']
            ic_type = 'isothermal_be' if 'isothermal' in suffix else 'lane_emden'

            ds_id = f"rp{orbit_info['r_peri']}pc"
            if suffix:
                ds_id += f"_{suffix}"
            ds_id += f"_{physics}"

            name_parts = [f"r_p={orbit_info['r_peri']}pc"]
            if 'isothermal' in suffix:
                name_parts.append("Isothermal BE")
            elif 'hires' in suffix:
                name_parts.append("Hi-Res")
            else:
                name_parts.append("Lane-Emden")
            name_parts.append("Cooling" if physics == 'cooling' else "Adiabatic")

            rel_path = os.path.relpath(results_dir, survey_dir / 'viz')

            imbh = config.get('imbh_parameters', {})
            cloud_pos = imbh.get('cloud_initial_position', {'0': 0, '1': 0, '2': 0})
            cloud_vel = imbh.get('cloud_initial_velocity', {'0': 0, '1': 0, '2': 0})

            dataset = {
                'id': ds_id,
                'name': ' '.join(name_parts),
                'path': rel_path,
                'orbit': orbit_info['orbit_id'],
                'physics': physics,
                'ic_type': ic_type,
                'n_snapshots': n_snapshots,
                'config': {
                    'r_peri': orbit_info['r_peri'],
                    'cloud_pos0': [cloud_pos.get('0', 0), cloud_pos.get('1', 0), cloud_pos.get('2', 0)],
                    'cloud_vel0': [cloud_vel.get('0', 0), cloud_vel.get('1', 0), cloud_vel.get('2', 0)],
                    'description': config.get('name', f'{" ".join(name_parts)} simulation')
                }
            }

            if ic_type == 'isothermal_be':
                for key in ['R_cloud', 'M_cloud', 'T_cloud', 'n_center']:
                    if key in config:
                        dataset['config'][key] = config[key]

            all_datasets.append(dataset)

    # Sort by orbit, ic_type, physics
    def sort_key(ds):
        orbit_num = float(ds['orbit'].replace('rp', '').replace('pc', ''))
        ic_priority = 0 if ds['ic_type'] == 'isothermal_be' else 1
        physics_priority = 0 if ds['physics'] == 'cooling' else 1
        return (orbit_num, ic_priority, physics_priority)

    all_datasets.sort(key=sort_key)

    return {
        'datasets': all_datasets,
        'common': {
            'G': 0.00430091,
            'M_BH': 100000.0,
            'tempConversion': 186.0,
            'densityToNH2': 20.3,
            'timeToMyr': 0.978
        },
        '_auto_discovery': True
    }


class AutoDiscoveryHandler(http.server.SimpleHTTPRequestHandler):
    """HTTP handler that auto-generates datasets.json on request."""

    def __init__(self, *args, survey_dir=None, **kwargs):
        self.survey_dir = survey_dir
        super().__init__(*args, **kwargs)

    def do_GET(self):
        # Auto-generate datasets.json when requested
        if self.path == '/viz/datasets.json' or self.path == '/datasets.json':
            self.send_datasets_json()
        else:
            super().do_GET()

    def send_datasets_json(self):
        """Generate and send datasets.json dynamically."""
        try:
            data = scan_datasets(self.survey_dir)
            content = json.dumps(data, indent=4).encode('utf-8')

            self.send_response(200)
            self.send_header('Content-Type', 'application/json')
            self.send_header('Content-Length', len(content))
            self.send_header('Cache-Control', 'no-cache')  # Always fresh
            self.end_headers()
            self.wfile.write(content)
        except Exception as e:
            self.send_error(500, f'Error generating datasets: {e}')


def main():
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8888

    # Determine survey directory
    script_dir = Path(__file__).parent.resolve()
    if script_dir.name == 'viz':
        survey_dir = script_dir.parent
    else:
        survey_dir = script_dir

    os.chdir(survey_dir)

    # Create handler with survey_dir
    handler = partial(AutoDiscoveryHandler, survey_dir=survey_dir)

    with socketserver.TCPServer(("", port), handler) as httpd:
        print(f"Auto-discovery server running at http://localhost:{port}/viz/")
        print(f"Scanning: {survey_dir}")
        print("datasets.json is generated dynamically on each request")
        print("Press Ctrl+C to stop")
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nServer stopped")


if __name__ == '__main__':
    main()
