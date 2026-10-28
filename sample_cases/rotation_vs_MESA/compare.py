"""Compare YREC and MESA runs of the solid-body star.

  python compare.py YREC_HISTORY MESA_HISTORY

YREC against MESA: a 1 Msun star rotating as a solid body, with no
wind and no mass loss, from the pre-main sequence to 1 Gyr, with the
same angular momentum J. Under J conservation omega = J/I, I = k^2 M R^2,
so the codes' rotation curves agree to the extent their I(t) does.
Prints each code's J conservation, then, against age, YREC's omega
over MESA's (as computed and scaled to the same J), I = J/omega, R, and
k^2 = J/(omega M R^2), with YREC's omega/omega_crit.
"""
import sys, numpy as np, pathlib
def load(p):
    lines = pathlib.Path(p).read_text().splitlines(); names = lines[5].split()
    a = np.array([list(map(float, l.split())) for l in lines[6:] if l.strip()])
    return {n: a[:, i] for i, n in enumerate(names)}
y = load(sys.argv[1]); m = load(sys.argv[2])
MSUN, RSUN = 1.989e33, 6.957e10
for name, h in (("YREC", y), ("MESA", m)):
    J = h["total_angular_momentum"]
    h["k2"] = J/(h["surf_avg_omega"]*MSUN*(10**h["log_R"]*RSUN)**2)
    print(f"{name}: {len(J)} models, log R from {h['log_R'][0]:.3f}, J = {J[0]:.4e}, "
          f"J/J0 - 1 in [{J.min()/J[0]-1:+.1e}, {J.max()/J[0]-1:+.1e}]")
print(f"J(YREC)/J(MESA) - 1 = {y['total_angular_momentum'][0]/m['total_angular_momentum'][0]-1:+.4f}")
mr = m["center_omega"]/m["surf_avg_omega"] - 1
print(f"MESA rigidity: |omega_c/omega_s - 1| <= {np.abs(mr).max():.1e}")
jr = y["total_angular_momentum"][0]/m["total_angular_momentum"][0]
print("\n   age (yr)  omega YREC  omega MESA  ratio-1  J-scaled   I ratio-1  R ratio-1  k2 YREC  k2 MESA  ratio-1  w/w_crit YREC")
G = 6.674e-8
for t in (3e6, 5e6, 1e7, 2e7, 3e7, 4e7, 5e7, 1e8, 3e8, 1e9):
    f = lambda h, c: np.interp(np.log10(t), np.log10(h["star_age"]), h[c])
    wy, wm = f(y, "surf_avg_omega"), f(m, "surf_avg_omega")
    ry, rm = 10**f(y, "log_R"), 10**f(m, "log_R")
    wc = np.sqrt(G*MSUN/(ry*RSUN)**3)
    iy = f(y, "total_angular_momentum")/wy; im = f(m, "total_angular_momentum")/wm
    print(f"  {t:9.2e}  {wy:10.4e}  {wm:10.4e}  {wy/wm-1:+.4f}  {wy/wm/jr-1:+.4f}    {iy/im-1:+.4f}    {ry/rm-1:+.4f}    {f(y,'k2'):7.4f}  {f(m,'k2'):7.4f}  {f(y,'k2')/f(m,'k2')-1:+.4f}   {wy/wc:.3f}")
