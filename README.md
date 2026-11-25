# Reciprocal sums of higher-order linear recurrences: effective rounding and the asymptotics of the error

Artatrana Suna and Prasanta Kumar Ray.

This repository contains the Python scripts, and the output they produced, for the numerical results of the paper
(Tables 1-5, the numerical observations (i)-(v) and the numerical remarks of the numerical-results section and of the section on small indices for the `k`-generalized Fibonacci numbers). All computations use **mpmath multiprecision
floating-point arithmetic**. No step uses interval arithmetic, so every numerical result is *computer-assisted*, not
formally certified (see "What is and is not certified" below). The theorems of the paper are proved for every Pisot
recurrence; the scripts verify finitely many cases numerically.

## Requirements

Python 3.9 or newer and `mpmath` 1.3 (the version used for the paper):

```bash
python -m pip install -r requirements.txt
```

## Scripts and what they reproduce

| Script | Reproduces | Output file | Approx. run time |
|--------|------------|-------------|------------------|
| `verify_pisot.py` | Table 1 (constants `alpha, mu, d_1, C_0, C_1`, index `N_0^root`, least starting index `n_*`), Table 2 (root-free index `N_0^exp` under Perron's condition), Table 3 (`alpha_2`, `kappa_2`, `N_par`) | `verify_pisot_output.txt`, `tables_generated.txt` | about 1 minute |
| `check_error_asymptotics.py` | Table 4 (exact `E_n` against the root sum), Table 5 (residual against the quadratic term), Tribonacci sign statistics and `theta/pi`, uniformity in `k` for `k`-generalized Fibonacci numbers, the closed form of `E_n` for `3 <= n <= k+1` and the value of `E_{k+2}` | `check_error_asymptotics_output.txt` | about 4 minutes |
| `verify_claims.py` | Remaining numerical claims of the numerical-results section: maxima of `abs(E_n)` for the eleven recurrences of Table 1, the floor rule for `2 <= n <= 59` in the seven second-order rows of Table 3, `abs(E_n) < 1/2` for 30-bonacci up to `n = 450`, the values of `(1-mu)k^3`, and the check that the non-dominant roots of maximal modulus form exactly one non-real conjugate pair for `F^(k)` (`4 <= k <= 40`), `P^(k)` (`4 <= k <= 15`) and `(4;1,0,0,1)` | `verify_claims_output.txt` | about 2 minutes |
| `verify_section6.py` | The values quoted at the end of the section on small indices for `F^(k)`: `(1/2 - x_k)4^k` for `k = 2, 3, 4, 6, 12, 30, 40`; `S_1^{-1}` for Tribonacci numbers (observation (i)); `theta/pi`, the count of positive `E_n` among `20 <= n <= 330` (157 of 311) and the maximal deviation of the normalized error from the leading cosine for `3 <= n <= 45` (observation (iv)); the constants `N_0^root (1-mu)` for `k = 8, 30` (remark on growth with `k`; `N_0^root = 14, 447` are taken from `verify_pisot_output.txt`) | `verify_section6_output.txt` | under 1 minute |

Run times are approximate. The scripts are deterministic: rerunning them with the same mpmath version gives the same
output as the files provided here.

## Usage

```bash
python verify_pisot.py                       # Tables 1-3; also writes tables_generated.txt
python verify_pisot.py --radii               # inclusion radii and margins 1-q for the eleven recurrences
python check_error_asymptotics.py            # error-expansion tables and observations
python verify_claims.py                      # remaining numerical claims; exit status 0 if all pass
python verify_claims.py --quick              # skips k = 160 and shortens the n <= 450 scan of 30-bonacci
python verify_section6.py                    # remaining values quoted in the paper; exit status 0 if all pass
```

To regenerate the output files:

```bash
python verify_pisot.py > verify_pisot_output.txt
python check_error_asymptotics.py > check_error_asymptotics_output.txt
python verify_claims.py > verify_claims_output.txt
python verify_section6.py > verify_section6_output.txt
```

The rows of Tables 1 and 2, in order, are: Fibonacci, Pell, (2;3,1), Tribonacci, (3;6,1,1), 8-bonacci, 30-bonacci,
(4;3,2,2,1), (5;2,2,1,1,1), Narayana (3;1,0,1), (4;1,0,0,1); the Perron-condition recurrences (3;3,0,1),
(4;3,0,0,1), (4;5,-1,0,1) are appended automatically.

## What is and is not certified

* **Controlled by precision margins and safety factors, not by interval arithmetic:** the roots of the characteristic
  polynomial (200 digits); the a-posteriori inclusion radii `k|Phi(z_i)|/prod|z_i - z_j|` (floating-point residuals; a
  printed radius of 0.0 means the residual underflowed at the working precision); the constants `C_0`, `C_1`, `q`, `d_1`
  (inflated or deflated by a relative `1e-30`); the thresholds `delta_n < 1`, `H_n < 1`, `B(n) < 1/2`; and the
  enclosures of `S_n` for `n < N_0` (working precision at least 60 digits beyond what is needed).
* **Not proofs:** the comparisons in `check_error_asymptotics.py` and `verify_claims.py`, which test the formulas of the
  paper at finitely many indices.

## Citation and licence

See `CITATION.cff`. Released under the MIT licence (see `LICENSE`).
