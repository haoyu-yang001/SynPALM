"""Completeness check of the revision simulations: expected vs present output files
(test arrays 999/9999 ignored) and unusually small files (possible truncated writes).
  python3 code/sh_run/check_sim_outputs.py
Prints, for each study, the number of missing outputs and the missing array IDs, so
that they can be resubmitted (all scripts skip outputs that already exist)."""
import os, statistics, collections
R = '<PROJECT_DIR>/simulation_code/results_rev1/'
def tl(d, names, nshow=15):
    have = set(os.listdir(R + d)) if os.path.isdir(R + d) else set()
    miss = sorted(n for n in names if n not in have)
    groups = collections.defaultdict(list)
    for n in names:
        if n in have: groups[n.split('_array_')[0]].append(os.path.getsize(R + d + n))
    small = []
    for n in names:
        if n in have:
            med = statistics.median(groups[n.split('_array_')[0]])
            if os.path.getsize(R + d + n) < 0.8 * med: small.append(n)
    print('%-22s expected %5d  present %5d  missing %4d  small %d' % (d, len(names), len(names) - len(miss), len(miss), len(small)))
    for n in miss[:nshow]: print('   missing', n)
    for n in small[:10]: print('   small  ', n)
t1 = [f'GRMr_{g}_{k}_missing_{m}_array_{a}.rds' for g, k in (('0.125', 'cousin'), ('0.25', 'halfsiblings'))
      for m in ('0', '0.25', '0.5', '0.75', '0.9') for a in range(1, 101)]
t1 += [f'GRMr_0.5_siblings_missing_{m}_array_{a}.rds' for m in ('0', '0.5', '0.75') for a in range(1, 501)]
t1 += [f'GRMr_0.5_siblings_missing_{m}_array_{a}.rds' for m in ('0.25', '0.9') for a in range(1, 101)]
tl('typeIerror/', t1)
tl('typeIerror_split/', [f'GRMr_0.5_split_missing_{m}_array_{a}.rds' for m in ('0.5', '0.75') for a in range(1, 201)])
tl('variance_components/', [f'missing_{m}_array_{a}.rds' for m in ('0', '0.25', '0.5', '0.75', '0.9') for a in range(1, 21)])
tl('S2/', [f'setting_{s:02d}_array_{a:03d}.rds' for s in range(1, 7) for a in range(1, 101)])
tl('S36/', [f'setting_{s:02d}_array_{a:03d}.rds' for s in range(1, 35) for a in range(1, 21)])
tl('S4/', [f'setting_{s:02d}_array_{a:03d}.rds' for s in range(1, 11) for a in range(1, 21)])
tl('S5/', [f'setting_{s:02d}_array_{a:03d}.rds' for s in range(1, 7) for a in range(1, 21)])
def rv(x):  # R's as.character() of a number, as in sprintf("%s", x)
    s = repr(float(x)).rstrip('0').rstrip('.') if float(x) != int(float(x)) else str(int(float(x)))
    return {'0.0001': '1e-04', '0.0002': '2e-04', '0.0003': '3e-04', '0.0004': '4e-04', '0.0005': '5e-04',
            '0.0006': '6e-04', '0.0007': '7e-04', '0.0008': '8e-04', '0.0009': '9e-04'}.get(s, s)
MS = ('0', '0.25', '0.5', '0.75', '0.9')
tl('true_geno_typeIerror_calgeno/', [f'rep_{k}_missing_{m}_rho_{r}_chunkarray_{a}_negative_log10_pvals_type1error_rev1.rds'
                                     for k in range(1, 18) for a in range(1, 41) for m in MS for r in ('0', '0.25', '0.5', '0.75')])
tl('true_geno_power_rho075/', [f'missing_{m}_rho_0.75_her_{rv(h / 1e4)}_chunkarray_{a}_negative_log10_pvals_for_power_rev1.rds'
                               for a in range(1, 1001) for m in MS for h in range(1, 8)])

## reruns of the settings kept from the submission (2026-10-16): power rho 0/0.25/0.5 (Fig 3, S5-S7),
## bias rho 0 (Fig 2A), true-genotype power rho 0/0.25/0.5 (Figs S11-S12)
KIN = (('0.125', 'cousin'), ('0.25', 'halfsiblings'), ('0.5', 'siblings'))
for r, tag in (('0', '000'), ('0.25', '025'), ('0.5', '050')):
    tl(f'power_rho{tag}/', [f'GRMr_{g}_missing_{m}_rho_{r}_rhogeno_{k}_her_{rv(h / 1e4)}_chunkarray_{a}_rev1.rds'
                            for g, k in KIN for m in MS for h in range(1, 11) for a in range(1, 201)], nshow=0)
tl('bias_rho000/', [f'GRMr_0.5_missing_{m}_rho_0_rhogeno_siblings_her_0_chunkarray_{a}_rev1.rds'
                    for m in MS for a in range(1, 1001)], nshow=0)
for r, tag in (('0', '000'), ('0.25', '025'), ('0.5', '050')):
    tl(f'true_geno_power_rho{tag}/', [f'missing_{m}_rho_{r}_her_{rv(h / 1e4)}_chunkarray_{a}_negative_log10_pvals_for_power_rev1.rds'
                                      for a in range(1, 1001) for m in MS for h in range(1, 8)], nshow=0)
