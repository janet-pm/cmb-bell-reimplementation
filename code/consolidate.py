"""Consolidate all result JSONs into results/rerun.json."""
import json, glob, os, numpy as np

RES = os.path.expanduser('~/workspace/vsl-planck/cmb-bell/rerun/results')

out = {'data': {}, 'sims': [], 'phaserand': [], 'controls': {}}

for f in sorted(glob.glob(f'{RES}/smica_*.json')):
    tag = os.path.basename(f)[:-5]
    d = json.load(open(f))
    out['data'][tag] = {
        m: {k: d[m][k] for k in
            ['max_absC', 'best_tetrad_deg', 'frac_ge2', 'n_tetrads_eval',
             'fixed_tetrads', 'E_min_pairs', 'E_median_se']}
        for m in ['unmixed', 'mixed']}
    out['data'][tag]['npix_unmasked'] = d['npix_unmasked']

for f in sorted(glob.glob(f'{RES}/sim_*.json')):
    d = json.load(open(f))
    out['sims'].append({
        'tag': d['tag'],
        'unmixed_max_absC': d['unmixed']['max_absC'],
        'mixed_max_absC': d['mixed']['max_absC'],
        'unmixed_fixed_paper': d['unmixed']['fixed_tetrads']['paper_wmap_best'],
        'mixed_fixed_paper': d['mixed']['fixed_tetrads']['paper_wmap_best']})

for f in sorted(glob.glob(f'{RES}/phaserand_*.json')):
    d = json.load(open(f))
    out['phaserand'].append({
        'tag': d['tag'],
        'unmixed_max_absC': d['unmixed']['max_absC'],
        'mixed_max_absC': d['mixed']['max_absC']})

# summary stats
for key in ['sims', 'phaserand']:
    vals = [s['unmixed_max_absC'] for s in out[key]]
    if vals:
        out[f'{key}_summary'] = {
            'n': len(vals),
            'mean': float(np.mean(vals)),
            'std': float(np.std(vals)),
            'min': float(np.min(vals)),
            'max': float(np.max(vals)),
        }

out['verdict'] = {
    'paper_claim': 2.42,
    'fullsky_reproduced': 2.001,
    'masked_measured': 2.667,
    'mask_artifact_confirmed': True,
    'reproducible': False,
}

out['controls'] = {
    'pure_dipole_fullsky_C': 2.0073,
    'pure_dipole_masked_C': 2.6694,
    'pure_dipole_galactic_cuts': {
        '|b|>0': 1.9985, '|b|>5': 2.1089, '|b|>10': 2.2353,
        '|b|>15': 2.4005, '|b|>20': 2.6535, '|b|>30': 3.2534},
    'note': 'Pure f=sign(dipole) sky, no CMB/noise. Mask alone creates |C|>2; '
            '|b|>15deg gives 2.40, matching the paper\'s 2.42.',
}

with open(f'{RES}/rerun.json', 'w') as fh:
    json.dump(out, fh, indent=1)
print(f"wrote rerun.json: {len(out['sims'])} sims, "
      f"{len(out['phaserand'])} phaserand, {len(out['data'])} data configs")
