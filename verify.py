#!/usr/bin/env python3
"""Independent verification script for the diagonal asymptotics of
OEIS A263160 (first-zero lattice paths) and OEIS A262810 (origin-reaching lattice paths).

Uses only the Python standard library. Performs:
  1. Direct combinatorial dynamic programming for small n vs. exact series formulas.
  2. High-precision evaluation of the closed-form constant C (50 significant digits).
  3. Comparison against official OEIS b-file values (optional CLI paths).
  4. Cancellation-free large-n numerical verification up to n = 1000.
"""

from __future__ import annotations

import argparse
import decimal
import math
from functools import lru_cache
from pathlib import Path


# ---------------------------------------------------------------------------
# 1. Direct combinatorial DP for small n (from literal definitions)
# ---------------------------------------------------------------------------

def exact_a_dp(n: int) -> int:
    """Exact A263160(n): paths from (n,...,n) stopping when any coordinate hits 0."""
    if n == 0:
        return 1

    @lru_cache(maxsize=None)
    def dp(state: tuple[int, ...]) -> int:
        total = 0
        for mask in range(1, 1 << n):
            nxt = []
            stopped = False
            for i in range(n):
                v = state[i] - ((mask >> i) & 1)
                if v == 0:
                    stopped = True
                    break
                nxt.append(v)
            if stopped:
                total += 1
            else:
                total += dp(tuple(sorted(nxt)))
        return total

    return dp((n,) * n)


def exact_b_dp(n: int) -> int:
    """Exact A262810(n): paths from (n,...,n) to (0,...,0) decrementing >=1 coord."""
    if n == 0:
        return 1

    @lru_cache(maxsize=None)
    def dp(state: tuple[int, ...]) -> int:
        if state[-1] == 0:
            return 1
        total = 0
        for mask in range(1, 1 << n):
            nxt = []
            valid = True
            for i in range(n):
                v = state[i] - ((mask >> i) & 1)
                if v < 0:
                    valid = False
                    break
                nxt.append(v)
            if valid:
                total += dp(tuple(sorted(nxt)))
        return total

    return dp((n,) * n)


# ---------------------------------------------------------------------------
# 2. Exact series evaluation for small n
# ---------------------------------------------------------------------------

def series_a_float(n: int) -> float:
    """Compute 1 + (2^(n-1)-1) * sum_{k>=0} 2^{-k} P_{n-1}(k)^n for small n."""
    if n <= 1:
        return 1.0
    m = n - 1
    s = 0.0
    for k in range(0, 1000):
        p_mk = sum(math.comb(k, r) for r in range(m + 1))
        term = (p_mk ** n) * (0.5 ** k)
        s += term
        if k > 20 and term < 1e-18 * s:
            break
    return 1.0 + (2 ** (n - 1) - 1) * s


def series_b_float(n: int) -> float:
    """Compute (1/2) * sum_{k>=n} 2^{-k} binom(k, n)^n for small n."""
    if n == 0:
        return 1.0
    s = 0.0
    for k in range(n, 1000):
        term = (math.comb(k, n) ** n) * (0.5 ** (k + 1))
        s += term
        if k > 20 and term < 1e-18 * s:
            break
    return s


# ---------------------------------------------------------------------------
# 3. High-precision closed-form constant C for A263160
# ---------------------------------------------------------------------------

def compute_C_decimal(digits: int = 50) -> decimal.Decimal:
    """Compute C = (2*sqrt(2*pi)/log(2)) * exp(5/12 - (log 2)^2 / 24)."""
    with decimal.localcontext() as ctx:
        ctx.prec = digits + 20

        def atan_inv(x: int) -> decimal.Decimal:
            inv_x = decimal.Decimal(1) / decimal.Decimal(x)
            inv_x2 = inv_x * inv_x
            term = inv_x
            total = term
            k = 1
            while True:
                term = -term * inv_x2
                add = term / decimal.Decimal(2 * k + 1)
                if abs(add) < decimal.Decimal(10) ** (-(digits + 18)):
                    break
                total += add
                k += 1
            return total

        pi = 16 * atan_inv(5) - 4 * atan_inv(239)
        ln2 = decimal.Decimal(2).ln()
        two = decimal.Decimal(2)
        exp_arg = decimal.Decimal(5) / decimal.Decimal(12) - (ln2 * ln2) / decimal.Decimal(24)
        c_val = (two * (two * pi).sqrt() / ln2) * exp_arg.exp()
    with decimal.localcontext() as ctx:
        ctx.prec = digits
        return +c_val


# ---------------------------------------------------------------------------
# 4. Cancellation-free large-n evaluation (avoids O(n^2 log n) float loss)
# ---------------------------------------------------------------------------

def _log1p_minus_u(u: float) -> float:
    """Accurately compute log(1 + u) - u for small |u|."""
    if abs(u) < 1e-2:
        return math.fsum(((-1) ** (d - 1)) * (u ** d) / d for d in range(2, 16))
    return math.log1p(u) - u


def _log_falling_ratio(k: int, m: int) -> float:
    """Compute sum_{j=1}^{m-1} log(1 - j/k) = log(binom(k, m) * m! / k^m)."""
    return math.fsum(math.log1p(-j / k) for j in range(1, m))


def _log_P_over_binom(k: int, m: int) -> float:
    """Compute log(P_m(k) / binom(k, m)) for k >= 2m."""
    ratio = 1.0
    term = 1.0
    for j in range(1, m + 1):
        term *= (m - j + 1) / (k - m + j)
        ratio += term
        if term < 1e-16 * ratio:
            break
    return math.log(ratio)


def eval_a_constant_stable(n: int, half_width_sigma: float = 12.0) -> float:
    """Evaluate C_n = a(n) * pi^(n/2) * (log 2)^(n(n-1)) / n^(n(n-3/2)+1)."""
    lam = math.log(2.0)
    m = n - 1
    N = n * m
    k0 = N / lam
    sigma = math.sqrt(N) / lam
    k_lo = max(2 * m, int(math.floor(k0 - half_width_sigma * sigma)))
    k_hi = int(math.ceil(k0 + half_width_sigma * sigma))

    r_m = math.lgamma(m + 1) - (0.5 * math.log(2.0 * math.pi * m) + m * math.log(m) - m)
    const_shift = (
        (0.5 * n - 1.0) * lam
        + math.log1p(-(0.5 ** (n - 1)))
        - math.log(n)
        - 0.5 * n * math.log1p(-1.0 / n)
        - n * r_m
    )

    log_terms = [
        N * _log1p_minus_u((k - k0) / k0)
        + n * (_log_falling_ratio(k, m) + _log_P_over_binom(k, m))
        for k in range(k_lo, k_hi + 1)
    ]
    max_lt = max(log_terms)
    total = math.fsum(math.exp(lt - max_lt) for lt in log_terms)
    return math.exp(const_shift + max_lt) * total


def eval_b_ratio_stable(n: int, half_width_sigma: float = 12.0) -> float:
    """Evaluate b(n) / b_asym(n) for OEIS A262810."""
    lam = math.log(2.0)
    m = n
    N = n * m
    k0 = N / lam
    sigma = math.sqrt(N) / lam
    k_lo = max(m, int(math.floor(k0 - half_width_sigma * sigma)))
    k_hi = int(math.ceil(k0 + half_width_sigma * sigma))

    r_n = math.lgamma(n + 1) - (0.5 * math.log(2.0 * math.pi * n) + n * math.log(n) - n)
    const_shift = (
        -math.log(n)
        + 1.0 / 12.0
        + (0.5 * n - 1.0 + lam / 24.0) * lam
        - 0.5 * math.log(math.pi)
        + math.log(lam)
        - n * r_n
    )

    log_terms = [
        N * _log1p_minus_u((k - k0) / k0) + n * _log_falling_ratio(k, m)
        for k in range(k_lo, k_hi + 1)
    ]
    max_lt = max(log_terms)
    total = math.fsum(math.exp(lt - max_lt) for lt in log_terms)
    return math.exp(const_shift + max_lt) * total


def log_a_scale(n: int) -> float:
    lam = math.log(2.0)
    return (
        (n * (n - 1.5) + 1.0) * math.log(n)
        - 0.5 * n * math.log(math.pi)
        - n * (n - 1) * math.log(lam)
    )


def log_b_asym(n: int) -> float:
    lam = math.log(2.0)
    return (
        (n * n - 0.5 * n + 1.0) * math.log(n)
        - 1.0 / 12.0
        - (n + lam / 24.0) * lam
        - 0.5 * (n - 1) * math.log(math.pi)
        - (n * n + 1) * math.log(lam)
    )


def load_bfile(path: Path) -> dict[int, int]:
    vals: dict[int, int] = {}
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        k_s, v_s = s.split()
        vals[int(k_s)] = int(v_s)
    return vals


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify A263160 and A262810 asymptotics.")
    parser.add_argument("--b263160", type=Path, default=None, help="Optional path to OEIS b263160.txt")
    parser.add_argument("--b262810", type=Path, default=None, help="Optional path to OEIS b262810.txt")
    args = parser.parse_args()

    print("=== 1. Exact Combinatorial DP vs. Geometric Binomial Series ===")
    for n in range(1, 5):
        a_dp = exact_a_dp(n)
        a_ser = series_a_float(n)
        b_dp = exact_b_dp(n)
        b_ser = series_b_float(n)
        assert abs(a_dp - a_ser) / max(1, a_dp) < 1e-12
        assert abs(b_dp - b_ser) / max(1, b_dp) < 1e-12
        print(f"n={n}: A263160={a_dp} (series={a_ser:.4f}), A262810={b_dp} (series={b_ser:.4f})")

    c_dec = compute_C_decimal(50)
    lam = math.log(2.0)
    c_flt = (2.0 * math.sqrt(2.0 * math.pi) / lam) * math.exp(5.0 / 12.0 - (lam * lam) / 24.0)
    print("\n=== 2. Closed-Form Constant C (50 significant digits) ===")
    print(f"C = {c_dec}")

    if args.b263160 and args.b263160.exists():
        print("\n=== 3a. Exact OEIS b263160.txt Comparison ===")
        bvals = load_bfile(args.b263160)
        for n in [5, 10, 15, 20, 25, 26]:
            if n in bvals:
                c_emp = math.exp(math.log(bvals[n]) - log_a_scale(n))
                rel = (c_emp - c_flt) / c_flt
                print(f"n={n:4d}: C_n = {c_emp:.10f}, rel_err = {rel*100:+.6f}%, n*rel_err = {n*rel:+.5f}")

    if args.b262810 and args.b262810.exists():
        print("\n=== 3b. Exact OEIS b262810.txt Comparison ===")
        bvals = load_bfile(args.b262810)
        for n in [5, 10, 15, 20, 25]:
            if n in bvals:
                ratio = math.exp(math.log(bvals[n]) - log_b_asym(n))
                rel = ratio - 1.0
                print(f"n={n:4d}: b(n)/asym = {ratio:.10f}, rel_err = {rel*100:+.6f}%, n^2*rel_err = {n*n*rel:+.5f}")

    print("\n=== 4a. Cancellation-Free Large-n Verification for A263160 ===")
    print(f"{'n':>6} | {'C_n (empirical)':>16} | {'Rel. error (%)':>14} | {'n * rel_err':>12}")
    for n in [10, 20, 26, 50, 100, 200, 500, 1000]:
        c_emp = eval_a_constant_stable(n)
        rel = (c_emp - c_flt) / c_flt
        print(f"{n:6d} | {c_emp:16.10f} | {rel*100:+14.6f} | {n*rel:+12.6f}")

    print("\n=== 4b. Cancellation-Free Large-n Verification for A262810 ===")
    print(f"{'n':>6} | {'b(n)/asym':>16} | {'Rel. error (%)':>14} | {'n^2 * rel_err':>14}")
    for n in [10, 20, 25, 50, 100, 200, 500, 1000]:
        ratio = eval_b_ratio_stable(n)
        rel = ratio - 1.0
        print(f"{n:6d} | {ratio:16.10f} | {rel*100:+14.6f} | {n*n*rel:+14.6f}")


if __name__ == "__main__":
    main()
