"""
verify_section6.py -- numerical statements of the paper that are not printed
by verify_pisot.py, check_error_asymptotics.py or verify_claims.py.

Reproduces
  (1) Section 6:  (1/2 - x_k) 4^k  for k = 2,3,4,6,12,30,40  (x_k = E_2 for F^(k))
  (2) Observation (i):  S_1^{-1} for Tribonacci numbers
  (3) Observation (iv): for Tribonacci, theta/pi, the number of n in [20,330]
      with E_n > 0 (out of 311), and the maximal deviation
      |E_n/(2|kappa|mu^(n-1)) - cos((n-1)theta + arg kappa)| over 3 <= n <= 45
  (4) Remark on growth: N_0^root * (1 - mu) for k = 8 and 30, with N_0^root
      taken from verify_pisot.py (14 and 447).

Requires only mpmath 1.3.  Run:  python verify_section6.py > verify_section6_output.txt
Exit status 0 if every check agrees with the value printed in the paper.
"""
import sys
from mpmath import mp, mpf, polyroots, cos, arg, nstr, pi

mp.dps = 120
ok = True


def check(label, computed, paper, tol):
    global ok
    good = abs(computed - paper) <= tol
    ok &= good
    print("[%s] %-44s computed %s, paper %s" %
          ("PASS" if good else "FAIL", label, nstr(computed, 6), paper))


def seq(c, N):
    """u_n for n = 2-k..N (dict), recurrence (k; c_1..c_k), u_1 = 1."""
    k = len(c)
    u = {i: mpf(0) for i in range(2 - k, 1)}
    u[1] = mpf(1)
    for n in range(2, N + 1):
        u[n] = sum(c[j - 1] * u[n - j] for j in range(1, k + 1))
    return u


def tail_inv(u, n, N):
    return 1 / sum(mpf(1) / u[m] for m in range(n, N + 1))


# (1) Section 6
print("1. Section 6: (1/2 - x_k) 4^k")
paper = {2: 1.220, 3: 0.9506, 4: 0.8428, 6: 0.7719, 12: 0.7503, 30: 0.75000, 40: 0.75000}
for k, p in paper.items():
    u = seq([1] * k, 400)
    x = tail_inv(u, 2, 400)          # E_2 = x_k since u_2 - u_1 = 0
    val = (mpf(1) / 2 - x) * 4 ** k
    check("(1/2 - x_k)4^k, k=%d" % k, val, p, mpf(5) * 10 ** (-(len(str(p)) - 2)) if p != 0.75 else mpf("5e-6"))

# (2)-(3) Tribonacci
print("\n2. Tribonacci, observation (i) and (iv)")
N = 800
u = seq([1, 1, 1], N)
check("S_1^{-1}, Tribonacci", tail_inv(u, 1, N), 0.3267, mpf("5e-5"))

rts = polyroots([1, -1, -1, -1], maxsteps=500, extraprec=400)
alpha = max(rts, key=abs)
a2 = [z for z in rts if abs(z - alpha) > mpf("1e-20") and z.imag > 0][0]
dj = a2 / (3 * a2 ** 2 - 2 * a2 - 1)               # alpha_j^(k-2)/Phi'(alpha_j), k=3
kap = dj * (alpha - a2) ** 2 / (alpha ** 2 - a2)
mu, th, w = abs(a2), arg(a2), arg(kap)
check("theta/pi", th / pi, 0.69272, mpf("5e-6"))

pos = sum(1 for n in range(20, 331) if tail_inv(u, n, N) - (u[n] - u[n - 1]) > 0)
print("   E_n > 0 for %d of %d indices n in [20,330] (proportion %s)" %
      (pos, 311, nstr(mpf(pos) / 311, 4)))
ok &= (pos == 157)

worst = mpf(0)
for n in range(3, 46):
    E = tail_inv(u, n, N) - (u[n] - u[n - 1])
    r = E / (2 * abs(kap) * mu ** (n - 1))
    worst = max(worst, abs(r - cos((n - 1) * th + w)))
print("   max deviation from leading cosine, 3<=n<=45: %s" % nstr(worst, 4))
ok &= (worst <= mpf("3e-3"))

# (4) growth remark
print("\n3. Remark on growth: N_0^root (1 - mu)")
for k, N0, p in ((8, 14, 0.74), (30, 447, 0.68)):
    r = polyroots([1] + [-1] * k, maxsteps=1000, extraprec=800)
    al = max(r, key=abs)
    mu_k = max(abs(z) for z in r if abs(z - al) > mpf("1e-20"))
    check("N0*(1-mu), k=%d (N0=%d)" % (k, N0), N0 * (1 - mu_k), p, mpf("5e-3"))

print("\n" + ("All checks passed." if ok else "SOME CHECKS FAILED."))
sys.exit(0 if ok else 1)
