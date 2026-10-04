"""Orchestrate the full clean-room rerun. Run steps individually; this driver
runs them in order and collects results/rerun.json.

Steps:
  1. prep.py            (needs PLA data files)
  2. analyze.py on SMICA neigh + harm      (headline: real data)
  3. analyze.py on HM-diff neigh + harm    (pure-noise diagnostic)
  4. sims.py --n 100                       (Gaussian null)
  5. sims.py --phase-rand --n 5            (phase-randomized diagnostic)
  6. collect -> results/rerun.json + summary print
"""
import json, os, subprocess, sys, glob

HERE = os.path.expanduser('~/workspace/vsl-planck/cmb-bell/rerun')
CODE = f'{HERE}/code'
DATA = f'{HERE}/data'
RES = f'{HERE}/results'
PY = f'{HERE}/venv/bin/python'


def run(cmd):
    print('+', ' '.join(cmd), flush=True)
    subprocess.run(cmd, check=True)


def main():
    os.makedirs(RES, exist_ok=True)
    run([PY, f'{CODE}/prep.py'])
    run([PY, f'{CODE}/analyze.py', f'{DATA}/prep_smica_neigh.npz',
         'smica_neigh', '100000000', '3000000', '0'])
    run([PY, f'{CODE}/analyze.py', f'{DATA}/prep_smica_harm.npz',
         'smica_harm', '100000000', '3000000', '1'])
    run([PY, f'{CODE}/analyze.py', f'{DATA}/prep_hmdiff_neigh.npz',
         'hmdiff_neigh', '30000000', '1500000', '2'])
    run([PY, f'{CODE}/analyze.py', f'{DATA}/prep_hmdiff_harm.npz',
         'hmdiff_harm', '30000000', '1500000', '3'])
    run([PY, f'{CODE}/sims.py', '--n', '100', '--n-pairs', '30000000',
         '--n-tetrads', '1500000'])
    run([PY, f'{CODE}/sims.py', '--phase-rand', '--n', '5', '--n-pairs',
         '30000000', '--n-tetrads', '1500000'])

    # collect
    rerun = {}
    for jf in sorted(glob.glob(f'{RES}/*.json')):
        if os.path.basename(jf) == 'rerun.json':
            continue
        with open(jf) as fh:
            rerun[os.path.basename(jf)[:-5]] = json.load(fh)
    with open(f'{RES}/rerun.json', 'w') as fh:
        json.dump(rerun, fh)
    print('=' * 60)
    for tag in ['smica_neigh', 'smica_harm', 'hmdiff_neigh', 'hmdiff_harm']:
        r = rerun[tag]
        for m in ['unmixed', 'mixed']:
            print(f"{tag}/{m}: max|C|={r[m]['max_absC']:.4f} "
                  f"fixed(paper)={r[m]['fixed_tetrads']['paper_wmap_best']:.4f} "
                  f"frac>=2={r[m]['frac_ge2']:.5f}")
    for h in ['neigh', 'harm']:
        mx = [rerun[f'sim_{i:03d}_{h}'][m]['max_absC']
              for i in range(100) for m in ['unmixed', 'mixed']]
        import numpy as np
        mx = np.array(mx)
        print(f"sims/{h}: max|C| median={np.median(mx):.4f} "
              f"p99={np.percentile(mx,99):.4f} max={mx.max():.4f} (n={len(mx)})")
    print('ALL DONE', flush=True)


if __name__ == '__main__':
    main()
