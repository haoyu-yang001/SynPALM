"""LaTeX rows of Table 1 (type I error of Syn-PALM and POP-GWAS, sib pairs) from
simulation_code/results_rev1/figs/table1_rev1.csv written by analyze_t1e_rev1.R.
  python3 simulation_code/rev1/table1_tex_rev1.py > simulation_code/results_rev1/figs/table1_rows_rev1.tex
Rows in the order of the submitted table: rho blocks 0, 0.25, 0.5, 0.75; within each, pi = 0, 0.5, 0.75."""
import csv
B = '<PROJECT_DIR>/simulation_code/results_rev1/figs/'
r = {(float(x['miss']), float(x['rho'])): x for x in csv.DictReader(open(B + 'table1_rev1.csv'))}
def f(v):
    v = float(v)
    return '0' if v == 0 else '%.1e' % v   # same style as the submitted table (1.0e-03)
out = []
for rho in (0, 0.25, 0.5, 0.75):
    for miss in (0, 0.5, 0.75):
        x = r[(miss, rho)]
        syn = [f(x['%s_SynSurrG' % a]) for a in ('1e-03', '1e-05', '1e-07')]
        pop = [f(x['%s_POPGWAS' % a]) for a in ('1e-03', '1e-05', '1e-07')]
        out.append('%g & %g & %s && %s \\\\' % (miss, rho, ' & '.join(syn), ' & '.join(pop)))
    out[-1] += '\\hline' if rho < 0.75 else ''
print('\n'.join(out))
