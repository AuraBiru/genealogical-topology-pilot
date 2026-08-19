#!/usr/bin/env python3
"""
Directional Coherence in Heterogeneous Cosmological Anomalies
==============================================================
Consolidated reproducible analysis.

Reproduces every statistic in the companion paper (Biru, Claude & Muse 2026b)
from data/verified_independent.csv. Seeded for stability; Monte Carlo values
match the paper within sampling noise.

Usage:  python code/analysis.py [--fast]
        --fast reduces Monte Carlo sizes for a quick check (~10x faster).

Outputs: results/results_summary.txt, results/sky_map.png,
         results/fanning_fit.png

Licence: CC-BY 4.0. Authors: Aura Biru, Claude (Anthropic), Muse (Meta).
"""

import csv
import sys
from itertools import combinations, permutations
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import rankdata, spearmanr

RNG = np.random.default_rng(7)
FAST = "--fast" in sys.argv
N_MC = 5_000 if FAST else 50_000
N_LE = 100 if FAST else 600
N_TRUNK = 2_000 if FAST else 8_000
ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
GA = (307.0, 9.0)
LOG_ZMAX = np.log10(1101)
LINES = []


def log(msg=""):
    print(msg)
    LINES.append(str(msg))


def angsep(l1, b1, l2, b2):
    l1, b1, l2, b2 = map(np.radians, [l1, b1, l2, b2])
    cs = np.sin(b1) * np.sin(b2) + np.cos(b1) * np.cos(b2) * np.cos(l1 - l2)
    return np.degrees(np.arccos(np.clip(cs, -1, 1)))


def load():
    core, contested = [], []
    with open(ROOT / "data" / "verified_independent.csv") as fh:
        for row in csv.DictReader(fh):
            entry = (row["name"], float(row["l_deg"]), float(row["b_deg"]), float(row["z"]))
            (core if row["status"] == "core" else contested).append(entry)
    return core, contested


def mean_pairwise_sep(pts):
    return float(np.mean([angsep(a[1], a[2], b[1], b[2]) for a, b in combinations(pts, 2)]))


def k_simple(a, b):
    th = angsep(a[1], a[2], b[1], b[2])
    dz = abs(np.log10(1 + a[3]) - np.log10(1 + b[3]))
    return 0.5 * th / 180 + 0.5 * dz / LOG_ZMAX


def random_sky(n):
    phi = RNG.uniform(0, 2 * np.pi, n)
    ct = RNG.uniform(-1, 1, n)
    return np.degrees(phi), 90 - np.degrees(np.arccos(ct))


def mc_pvalue(pts, statistic, observed, n_iter):
    """One-sided p for statistic <= observed under random sky positions."""
    zs = [p[3] for p in pts]
    count = 0
    for _ in range(n_iter):
        rl, rb = random_sky(len(pts))
        rand_pts = [("r", rl[i], rb[i], zs[i]) for i in range(len(pts))]
        if statistic(rand_pts) <= observed:
            count += 1
    return count / n_iter


def exact_perm_p(pts, ref=GA):
    dev = np.array([angsep(p[1], p[2], *ref) for p in pts])
    lz = np.array([np.log10(1 + p[3]) for p in pts])
    zr = rankdata(lz)
    r_obs = spearmanr(lz, dev)[0]
    count = total = 0
    for perm in permutations(range(len(pts))):
        r = spearmanr(zr, dev[list(perm)])[0]
        if r >= r_obs - 1e-12:
            count += 1
        total += 1
    return r_obs, count / total


def unit(l, b):
    l, b = np.radians(l), np.radians(b)
    return np.stack([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)], -1)


def sky_scan(pts, grid_step=3):
    gl = np.arange(0, 360, grid_step)
    gb = np.arange(-87, 90, grid_step)
    GL, GB = np.meshgrid(gl, gb)
    grid = unit(GL.ravel(), GB.ravel())
    P = unit(np.array([p[1] for p in pts]), np.array([p[2] for p in pts]))
    lz = np.array([np.log10(1 + p[3]) for p in pts])
    zr = rankdata(lz)
    zc = zr - zr.mean()
    zn = np.sqrt((zc ** 2).sum())

    def max_rho(points):
        dev = np.degrees(np.arccos(np.clip(grid @ points.T, -1, 1)))
        rk = np.argsort(np.argsort(dev, axis=1), axis=1).astype(float)
        rc = rk - rk.mean(axis=1, keepdims=True)
        rho = (rc @ zc) / (np.sqrt((rc ** 2).sum(axis=1)) * zn)
        i = int(np.argmax(rho))
        return float(rho[i]), i

    r_obs, i_obs = max_rho(P)
    exceed = 0
    for _ in range(N_LE):
        rl, rb = random_sky(len(pts))
        r, _ = max_rho(unit(rl, rb))
        if r >= r_obs:
            exceed += 1
    return r_obs, (float(GL.ravel()[i_obs]), float(GB.ravel()[i_obs])), exceed / N_LE


def f_log(d, a, b):
    return a * np.log(1 + d) + b


def f_lin(d, a, b):
    return a * d + b


def f_sig(d, L, k, d0):
    return L / (1 + np.exp(-k * (d - d0)))


def aic(y, yhat, npar):
    n = len(y)
    return n * np.log(np.sum((y - yhat) ** 2) / n) + 2 * npar


def four_point(pts, tol=0.05):
    n = len(pts)
    K = np.zeros((n, n))
    for i, j in combinations(range(n), 2):
        K[i][j] = K[j][i] = k_simple(pts[i], pts[j])
    passes = total = 0
    for A, B, C, D in combinations(range(n), 4):
        sums = sorted([K[A][B] + K[C][D], K[A][C] + K[B][D], K[A][D] + K[B][C]])
        if sums[2] > 0 and (sums[2] - sums[1]) / sums[2] <= tol:
            passes += 1
        total += 1
    return passes, total


def four_point_mc(pts, observed_passes, n_iter):
    zs = [p[3] for p in pts]
    exceed = 0
    rates = []
    for _ in range(n_iter):
        rl, rb = random_sky(len(pts))
        rand_pts = [("r", rl[i], rb[i], zs[i]) for i in range(len(pts))]
        p, _ = four_point(rand_pts)
        rates.append(p)
        if p >= observed_passes:
            exceed += 1
    return float(np.mean(rates)), exceed / n_iter


def main():
    core, contested = load()
    with_grb = core + contested

    log("=" * 64)
    log("VERIFIED CATALOGUE ANALYSIS  (seed=7)")
    log(f"Monte Carlo sizes: clustering/K {N_MC}, look-elsewhere {N_LE}")
    log("=" * 64)

    # 1. Clustering
    ms = mean_pairwise_sep(core)
    p_clu = mc_pvalue(core, mean_pairwise_sep, ms, N_MC)
    log(f"\n[1] Clustering (6-pt core): mean sep {ms:.1f} deg, p = {p_clu:.4f}")

    # 2. K_simple
    def mean_k(pts):
        return float(np.mean([k_simple(a, b) for a, b in combinations(pts, 2)]))

    mk = mean_k(core)
    p_k = mc_pvalue(core, mean_k, mk, N_MC)
    log(f"[2] K_simple (6-pt core): mean {mk:.3f}, p = {p_k:.4f}")

    # 3. Fanning, exact permutation
    r6, p6 = exact_perm_p(core)
    r7, p7 = exact_perm_p(with_grb)
    dev_core = [round(angsep(p[1], p[2], *GA)) for p in sorted(core, key=lambda x: x[3])]
    log(f"[3] Fanning from GA: core rho = {r6:.3f}, exact p = {p6:.4f}")
    log(f"    with contested GRB: rho = {r7:.3f}, exact p = {p7:.4f}")
    log(f"    deviation sequence (z-ordered, core): {dev_core}")

    # 4. Sky scan + look-elsewhere
    r_opt, opt, p_le = sky_scan(core)
    log(f"[4] Sky scan optimum: (l={opt[0]:.0f}, b={opt[1]:.0f}), rho = {r_opt:.3f}, "
        f"look-elsewhere p = {p_le:.2f}")

    # 5. Fits, AIC, z=8 predictions
    dev = np.array([angsep(p[1], p[2], *GA) for p in core])
    lz = np.array([np.log10(1 + p[3]) for p in core])
    po_log, pc_log = curve_fit(f_log, lz, dev, p0=[40, 20], maxfev=40000)
    e_log = np.sqrt(np.diag(pc_log))
    po_lin, _ = curve_fit(f_lin, lz, dev, p0=[30, 10], maxfev=40000)
    po_sig, _ = curve_fit(f_sig, lz, dev, p0=[100, 3, 0.3], maxfev=80000)
    a_l = aic(dev, f_log(lz, *po_log), 2)
    a_n = aic(dev, f_lin(lz, *po_lin), 2)
    a_s = aic(dev, f_sig(lz, *po_sig), 3)
    log(f"[5] Log fit: a = {po_log[0]:.1f} ± {e_log[0]:.1f}, b = {po_log[1]:.1f} ± {e_log[1]:.1f}")
    log(f"    Sigmoid: L = {po_sig[0]:.1f}, k = {po_sig[1]:.1f}, d0 = {po_sig[2]:.2f}")
    log(f"    AIC: sigmoid {a_s:.1f}, log {a_l:.1f}, linear {a_n:.1f} "
        f"(dAIC sig-log = {a_s - a_l:.1f})")
    d8 = np.log10(9)
    log(f"    z=8 predictions from GA: sigmoid {f_sig(d8, *po_sig):.0f}, "
        f"log {f_log(d8, *po_log):.0f}, linear {f_lin(d8, *po_lin):.0f} deg")

    # 6. Trunk with full covariance propagation
    a, b = po_log
    d180 = np.exp((180 - b) / a) - 1
    ef_central = (d180 - LOG_ZMAX) / np.log10(np.e)
    draws = RNG.multivariate_normal(po_log, pc_log, N_TRUNK)
    efs = []
    for aa, bb in draws:
        if aa <= 0:
            continue
        dd = np.exp((180 - bb) / aa) - 1
        efs.append((dd - LOG_ZMAX) / np.log10(np.e))
    efs = np.array(efs)
    log(f"[6] Trunk (log model): {ef_central:.0f} e-folds, "
        f"16-84% [{np.percentile(efs, 16):.0f}, {np.percentile(efs, 84):.0f}]")

    # 7. Four-point condition, verified sets
    p6q, t6 = four_point(core)
    rate6, pv6 = four_point_mc(core, p6q, N_MC // 10)
    p7q, t7 = four_point(with_grb)
    rate7, pv7 = four_point_mc(with_grb, p7q, N_MC // 10)
    pv6s = f"{pv6:.4f}" if pv6 > 0 else f"<{1/(N_MC//10):.4f}"
    pv7s = f"{pv7:.4f}" if pv7 > 0 else f"<{1/(N_MC//10):.4f}"
    log(f"[7] Four-point on K_simple: core {p6q}/{t6} pass "
        f"(random mean {rate6:.1f}), p = {pv6s}")
    log(f"    with contested GRB: {p7q}/{t7} pass (random mean {rate7:.1f}), p = {pv7s}")

    # Figures
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(9, 4.5))
        ax = fig.add_subplot(111, projection="mollweide")
        for name, l, bb_, z in with_grb:
            lr = np.radians(((l + 180) % 360) - 180) * -1
            br = np.radians(bb_)
            style = dict(s=70, edgecolors="k", zorder=3)
            colour = "tab:red" if z >= 1000 else ("tab:orange" if z >= 0.5 else "tab:blue")
            marker = "^" if name.startswith("GRB") else "o"
            ax.scatter(lr, br, c=colour, marker=marker, **style)
            ax.annotate(name.split()[0], (lr, br), fontsize=7,
                        xytext=(4, 4), textcoords="offset points")
        gl = np.radians(((GA[0] + 180) % 360) - 180) * -1
        ax.scatter(gl, np.radians(GA[1]), c="none", s=260, edgecolors="green",
                   linewidths=2, zorder=2, label="GA reference")
        ax.grid(alpha=0.3)
        ax.set_title("Verified anomaly catalogue (blue z<0.5, orange mid-z, red CMB; triangle = contested)")
        fig.tight_layout()
        fig.savefig(OUT / "sky_map.png", dpi=160)

        fig2, ax2 = plt.subplots(figsize=(7, 4.5))
        xs = np.linspace(0.0, 3.1, 300)
        ax2.plot(xs, f_sig(xs, *po_sig), label="sigmoid (AIC-preferred, core)")
        ax2.plot(xs, f_log(xs, *po_log), label="logarithmic")
        ax2.plot(xs, f_lin(xs, *po_lin), label="linear")
        ax2.scatter(lz, dev, c="k", zorder=3, label="verified core")
        g = contested[0]
        ax2.scatter([np.log10(1 + g[3])], [angsep(g[1], g[2], *GA)],
                    marker="^", c="tab:orange", edgecolors="k", zorder=3,
                    label="GRB (contested)")
        ax2.set_xlabel("log10(1+z)")
        ax2.set_ylabel("deviation from GA (deg)")
        ax2.legend(fontsize=8)
        ax2.set_title("Fanning: angular deviation vs depth")
        fig2.tight_layout()
        fig2.savefig(OUT / "fanning_fit.png", dpi=160)
        log("\nFigures written to results/.")
    except Exception as exc:  # matplotlib optional
        log(f"\nFigure generation skipped: {exc}")

    (OUT / "results_summary.txt").write_text("\n".join(LINES) + "\n")
    log("\nSummary written to results/results_summary.txt")


if __name__ == "__main__":
    main()
