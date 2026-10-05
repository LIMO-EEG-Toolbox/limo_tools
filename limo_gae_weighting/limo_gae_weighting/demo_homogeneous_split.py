#!/usr/bin/env python3
"""One homogeneous cohort -> random split (true difference = 0) -> known manipulations ->
compare LIMO Mean / Trimmed Mean (Yuen 20%) / weighted. Prints a JSON report.

    python demo_homogeneous_split.py                                        # inert?
    python demo_homogeneous_split.py --effect-amp 2.5                       # retained?
    python demo_homogeneous_split.py --contam-frac 0.3 --contam-kind offset # fabricated?
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from simulate_group_betas import preset, simulate, effect_map          # noqa: E402
from gae_weighting import (compute_weights, weighted_welch_t, yuen_t,   # noqa: E402
                           build_knn_adjacency)

BETA = 0


def inject_effect(datain, arm1, amp, coords, times):
    cfg = preset("group_effect", effect_amp=amp)
    eff = effect_map(cfg, coords, times)
    mask = np.abs(eff) >= cfg.effect_mask_frac * np.abs(eff).max()
    d = datain.copy()
    d[BETA][:, :, arm1] += eff.T[:, :, None]
    return d, eff, mask


def inject_contamination(datain, arm1, frac, kind, coords, rng, k_chan=6, amp=15.0):
    """offset = shared frontal deposit (fabricates); noise = heavy channel noise (swamps)."""
    d = datain.copy()
    nChan, nSub = d.shape[2], d.shape[3]
    carriers = rng.choice(arm1, max(1, int(round(frac * len(arm1)))), replace=False)
    truth = np.zeros((nChan, nSub), bool)
    frontal = np.argsort(-coords[:, 1])[:k_chan]
    for s in carriers:
        if kind == "offset":
            ch = frontal
            d[:, :, ch, s] += amp
        else:
            ch = rng.choice(nChan, k_chan, replace=False)
            d[:, :, ch, s] += rng.normal(0, amp, d[:, :, ch, s].shape)
        truth[ch, s] = True
    return d, truth, carriers


def perm_fwer(stat_fn, Y, groups, n_perm, rng, observed_max):
    mx = np.array([np.abs(stat_fn(Y, rng.permutation(groups))).max() for _ in range(n_perm)])
    return float((np.sum(mx >= observed_max) + 1) / (n_perm + 1)), float(np.quantile(mx, 0.95))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--n-per-group", type=int, default=20)
    p.add_argument("--effect-amp", type=float, default=0.0)
    p.add_argument("--contam-frac", type=float, default=0.0)
    p.add_argument("--contam-kind", default="offset", choices=["offset", "noise"])
    p.add_argument("--contam-amp", type=float, default=15.0)
    p.add_argument("--tau", type=float, default=3.0)
    p.add_argument("--seeds", type=int, default=1)
    p.add_argument("--n-perm", type=int, default=200)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--knn", type=int, default=6)
    p.add_argument("--json", default=None)
    a = p.parse_args()
    rng = np.random.default_rng(a.seed)

    out = simulate(preset("null", n_per_group=2 * a.n_per_group, n_groups=1, seed=a.seed))
    datain, coords, times = out["datain"], out["coords"], out["times"]
    nSub = datain.shape[3]
    groups = np.zeros(nSub, int)
    groups[rng.permutation(nSub)[:a.n_per_group]] = 1
    arm1 = np.where(groups == 1)[0]
    A = build_knn_adjacency(coords, a.knn) if a.knn > 0 else out["adjacency"]

    eff_mask = np.zeros((datain.shape[2], datain.shape[1]), bool)
    eff = None
    if a.effect_amp > 0:
        datain, eff, eff_mask = inject_effect(datain, arm1, a.effect_amp, coords, times)
    truth = np.zeros((datain.shape[2], nSub), bool)
    carriers = np.array([], int)
    if a.contam_frac > 0:
        datain, truth, carriers = inject_contamination(datain, arm1, a.contam_frac,
                                                       a.contam_kind, coords, rng, amp=a.contam_amp)

    W = compute_weights(datain, A, beta=BETA, tau=a.tau, seeds=tuple(range(a.seeds)))
    Y = datain[BETA]
    W1 = np.ones_like(W)
    crit = stats.t.ppf(0.975, df=nSub - 2)

    est = {"mean":   lambda Yy, gg: weighted_welch_t(Yy, gg, W1)[0],
           "robust": lambda Yy, gg: yuen_t(Yy, gg)[0],
           "ours":   lambda Yy, gg: weighted_welch_t(Yy, gg, W)[0]}
    diff = {"mean": weighted_welch_t(Y, groups, W1)[2],
            "robust": yuen_t(Y, groups)[1],
            "ours": weighted_welch_t(Y, groups, W)[2]}

    flagged = W < 1
    rep = dict(
        setting=vars(a), n_sub=int(nSub), n_carriers=int(len(carriers)),
        weights=dict(
            flag_rate_arm0=float(flagged[:, groups == 0].mean()),
            flag_rate_arm1=float(flagged[:, groups == 1].mean()),
            flag_rate_arms_p=float(stats.mannwhitneyu(
                flagged[:, groups == 0].mean(0), flagged[:, groups == 1].mean(0)).pvalue),
            min_weight=float(W.min()),
            hit_rate_on_corrupted=(float(flagged[truth].mean()) if truth.any() else None),
            false_flag_rate_on_clean=float(flagged[~truth].mean())),
        estimators={})
    for k, fn in est.items():
        t = fn(Y, groups)
        r = dict(n_sig=int((np.abs(t) > crit).sum()),
                 n_sig_outside_effect=int(((np.abs(t) > crit) & ~eff_mask).sum()),
                 max_abs_t=float(np.abs(t).max()))
        if eff is not None:
            r["effect_retained"] = float(diff[k][eff_mask].mean() / eff[eff_mask].mean())
        pv, thr = perm_fwer(fn, Y, groups, a.n_perm, np.random.default_rng(a.seed + 7),
                            r["max_abs_t"])
        r.update(perm_fwer_p=pv, n_sig_fwer=int((np.abs(t) > thr).sum()))
        rep["estimators"][k] = r

    print(json.dumps(rep, indent=2))
    if a.json:
        with open(a.json, "w") as f:
            json.dump(rep, f, indent=2)


if __name__ == "__main__":
    main()
