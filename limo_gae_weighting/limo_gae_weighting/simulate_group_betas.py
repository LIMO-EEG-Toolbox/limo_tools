#!/usr/bin/env python3
"""
simulate_group_betas.py -- synthetic group-structured beta-map generator for
validating the group-conditional GAE weighting (the "simulate groups with a
measurable difference" testbed Cyril Pernet asked for).

WHY THIS EXISTS
---------------
ERP CORE and ds002718 are single-cohort healthy-adult datasets: there is NO
natural patient-vs-control grouping to validate a *between-group* weighting
method against. So we manufacture groups whose difference is KNOWN, and whose
contamination is KNOWN, and check that a quality-weighting method:
  (1) PRESERVES a real group effect (does not down-weight it as "outlier"),
  (2) down-weights ONLY genuinely corrupted subjects/channels,
  (3) keeps Type-I error controlled under the null.

This maps 1:1 onto the validation scenarios in:
  GROUP_AE_BEST_PRACTICE.md   Section 6
  GROUP_AE_BRAINSTORM.md      Section 5
  GROUP_AE_DESIGN_DIALOGUE.md Section 1.2.4 / D5

OUTPUT FORMAT (drop-in for the existing pipeline)
-------------------------------------------------
The core tensor `datain` has shape [nBeta, nTime, nChan, nSubject] -- exactly
what limo_group_outliers.m assembles and NiPyAEoutliers.py consumes. Adjacency
is [nChan, nChan] binary symmetric. We also save the GROUND TRUTH that only a
simulation can give you: the injected effect mask, which subjects/channels were
contaminated, the group labels, and the per-subject measurement error (SME).

Design decisions worth knowing:
  * beta 0 carries the injected group effect; beta 1 is a matched control
    condition with NO group effect (a built-in specificity check).
  * "individual differences" (amplitude jitter, latency shift, smooth subject
    fields) are SIGNAL and must survive weighting -- they are generated
    separately from "contamination", which is the only thing a good weighter
    should suppress. Latency shift in particular is signal, not noise (RIDE).
  * measurement noise sd is the per-subject SME (SD_single_trial / sqrt(n_trials)),
    so the precision-weight / inverse-variance idea in the design docs has a
    concrete anchor here.

Author: generated for Renxiang Qiu's LIMO x GAE project.
Runs under the `limo-gae` conda env (numpy only for generation).
"""

from __future__ import annotations

import argparse
import json
import os
from dataclasses import dataclass, field, asdict

import numpy as np

# --------------------------------------------------------------------------- #
# 1) Montage: a ~32-channel 10-20-style layout on a normalized 2D scalp disc.
#    x = right(+)/left(-), y = front(+)/back(-). Approximate but topologically
#    correct, which is all a controllable spatial testbed needs. Occipito-
#    temporal sites (PO7/PO8/O1/O2) are present so an N170-like posterior effect
#    can be placed at an anatomically sensible location.
# --------------------------------------------------------------------------- #
_MONTAGE_2D = {
    "Fp1": (-0.30, 0.95), "Fp2": (0.30, 0.95),
    "F7": (-0.80, 0.60), "F3": (-0.40, 0.60), "Fz": (0.00, 0.62),
    "F4": (0.40, 0.60), "F8": (0.80, 0.60),
    "FC5": (-0.60, 0.35), "FC1": (-0.20, 0.35), "FC2": (0.20, 0.35), "FC6": (0.60, 0.35),
    "T7": (-0.95, 0.00), "C3": (-0.45, 0.00), "Cz": (0.00, 0.00),
    "C4": (0.45, 0.00), "T8": (0.95, 0.00),
    "CP5": (-0.60, -0.35), "CP1": (-0.20, -0.35), "CP2": (0.20, -0.35), "CP6": (0.60, -0.35),
    "P7": (-0.80, -0.60), "P3": (-0.40, -0.60), "Pz": (0.00, -0.62),
    "P4": (0.40, -0.60), "P8": (0.80, -0.60),
    "PO7": (-0.55, -0.80), "PO3": (-0.25, -0.82), "PO4": (0.25, -0.82), "PO8": (0.55, -0.80),
    "O1": (-0.30, -0.95), "Oz": (0.00, -0.98), "O2": (0.30, -0.95),
}


def default_montage():
    names = list(_MONTAGE_2D.keys())
    coords = np.array([_MONTAGE_2D[n] for n in names], dtype=np.float64)
    return names, coords


def build_adjacency(coords: np.ndarray, radius: float = 0.42) -> np.ndarray:
    """Binary symmetric neighbour matrix (no self-loops) from 2D distances.

    `radius` is tuned so mean degree lands ~4-6, similar to a real EEG
    neighbouring matrix. Matches the [nChan,nChan] adjacency the pipeline uses.
    """
    n = coords.shape[0]
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    A = (d <= radius).astype(np.float32)
    np.fill_diagonal(A, 0.0)
    return A


# --------------------------------------------------------------------------- #
# 2) ERP building blocks: separable spatiotemporal components.
# --------------------------------------------------------------------------- #
def temporal_gaussian(times: np.ndarray, peak_ms: float, width_ms: float) -> np.ndarray:
    return np.exp(-0.5 * ((times - peak_ms) / width_ms) ** 2)


def spatial_bump(coords: np.ndarray, focus_xy, spread: float) -> np.ndarray:
    focus = np.asarray(focus_xy, dtype=np.float64)
    d2 = np.sum((coords - focus[None, :]) ** 2, axis=1)
    return np.exp(-0.5 * d2 / spread ** 2)


@dataclass
class Component:
    """One separable ERP component: amp * spatial_bump ⊗ temporal_gaussian."""
    name: str
    focus: tuple
    spread: float
    peak_ms: float
    width_ms: float
    amp: float  # microvolt-ish; sign encodes polarity

    def render(self, coords, times, latency_shift_ms=0.0):
        s = spatial_bump(coords, self.focus, self.spread)          # [nChan]
        t = temporal_gaussian(times, self.peak_ms + latency_shift_ms, self.width_ms)  # [nTime]
        return self.amp * np.outer(s, t)                            # [nChan, nTime]


# Canonical component library (rough visual-ERP morphology).
def _component_library():
    return {
        "P1":   Component("P1",  focus=(0.0, -0.90), spread=0.45, peak_ms=110, width_ms=35, amp=+2.5),
        "N170l": Component("N170l", focus=(-0.55, -0.80), spread=0.35, peak_ms=170, width_ms=30, amp=-3.5),
        "N170r": Component("N170r", focus=(0.55, -0.80), spread=0.35, peak_ms=170, width_ms=30, amp=-3.5),
        "P3":   Component("P3",  focus=(0.0, -0.35), spread=0.55, peak_ms=380, width_ms=90, amp=+4.0),
        "N2f":  Component("N2f", focus=(0.0, 0.55), spread=0.50, peak_ms=250, width_ms=50, amp=-2.0),
    }


def base_template(coords, times, nBeta: int) -> np.ndarray:
    """Grand-average template per beta condition: [nBeta, nChan, nTime].

    beta 0 = a face-like ERP (P1 + bilateral N170 + P3).
    beta 1 = a different condition (P1 + frontal N2 + smaller P3) with NO N170,
             so it serves as a matched control that should stay null under a
             posterior N170-focused group manipulation.
    Further betas recycle beta 0's morphology (deterministic).
    """
    lib = _component_library()
    templates = []
    recipes = [
        ["P1", "N170l", "N170r", "P3"],       # beta 0
        ["P1", "N2f", "P3"],                  # beta 1 (control condition)
    ]
    for b in range(nBeta):
        recipe = recipes[b] if b < len(recipes) else recipes[0]
        Bmap = np.zeros((coords.shape[0], times.shape[0]), dtype=np.float64)
        for cname in recipe:
            Bmap += lib[cname].render(coords, times)
        templates.append(Bmap)
    return np.stack(templates, axis=0)


# --------------------------------------------------------------------------- #
# 3) Scenario configuration.
# --------------------------------------------------------------------------- #
@dataclass
class ScenarioConfig:
    scenario: str = "group_effect"
    n_per_group: int = 20
    n_groups: int = 2
    nBeta: int = 2
    t_start_ms: float = -200.0
    t_end_ms: float = 600.0
    nTime: int = 120

    # normal individual differences (SIGNAL -- must survive weighting)
    amp_cv: float = 0.15          # component-amplitude coefficient of variation
    latency_sd_ms: float = 8.0    # per-subject latency jitter
    indiv_field_amp: float = 0.6  # smooth per-subject spatial/temporal field
    indiv_field_rank: int = 3

    # measurement noise -> SME anchor
    n_trials_mean: int = 120
    n_trials_sd: int = 25
    single_trial_sd: float = 12.0  # single-trial uV sd -> SME = sd/sqrt(n_trials)
    noise_spatial_smooth: float = 0.0  # 0 = white; >0 = spatially correlated

    # injected group effect (beta 0 only)
    effect_focus: tuple = (0.55, -0.80)  # near PO8 (posterior-lateral, N170-like)
    effect_focus2: tuple = (-0.55, -0.80)  # near PO7 (bilateral)
    effect_spread: float = 0.33
    effect_peak_ms: float = 200.0
    effect_width_ms: float = 35.0
    effect_amp: float = 0.0        # uV added to group>=1 mean at the focus
    effect_mask_frac: float = 0.30  # mask = |effect| > frac * max|effect|

    # contamination (QUALITY problems -- SHOULD be down-weighted)
    contam_type: str = "none"      # none|bad_channel|time_burst|variance_inflation
    contam_frac: float = 0.30      # fraction of subjects contaminated
    contam_which_group: int = -1   # -1 = any group; else restrict to that group
    bad_channel_k: int = 3         # channels affected per contaminated subject
    bad_channel_gain: float = 6.0  # noise-sd multiplier on bad channels
    burst_peak_ms: float = 300.0
    burst_width_ms: float = 40.0
    burst_amp: float = 15.0
    var_inflation_gain: float = 4.0

    # group-level noise (one whole group noisier -- within-group calibration test)
    group_noise_gain: float = 1.0  # applied to group 1's SME (1.0 = off)

    seed: int = 0


def preset(scenario: str, **overrides) -> ScenarioConfig:
    """Named presets that map to the validation protocol in the design docs."""
    base = dict(scenario=scenario)
    if scenario == "null":
        # two groups, identical distribution -> Type-I / FWER check
        base.update(effect_amp=0.0, contam_type="none", group_noise_gain=1.0)
    elif scenario == "group_effect":
        # a real, preserved between-group effect at a posterior cluster
        base.update(effect_amp=2.5, contam_type="none")
    elif scenario == "bad_channel":
        # localized quality problem in a subset of subjects
        base.update(effect_amp=2.5, contam_type="bad_channel", contam_frac=0.30)
    elif scenario == "time_burst":
        base.update(effect_amp=2.5, contam_type="time_burst", contam_frac=0.30)
    elif scenario == "variance_inflation":
        # first-level estimation noisier for a subset (higher SME)
        base.update(effect_amp=2.5, contam_type="variance_inflation", contam_frac=0.30)
    elif scenario == "latency_shift":
        # BIOLOGICAL individual difference (signal, not noise). One group is
        # latency-shifted; a good weighter must NOT down-weight it.
        base.update(effect_amp=0.0, contam_type="none", latency_sd_ms=25.0)
    elif scenario == "group_noise":
        # one whole group noisier overall, same mean -> within-group calibration
        base.update(effect_amp=2.5, contam_type="none", group_noise_gain=2.5)
    elif scenario == "mixed":
        # realistic: real effect + contamination in a subset of one group
        base.update(effect_amp=2.5, contam_type="bad_channel", contam_frac=0.25,
                    contam_which_group=1)
    else:
        raise ValueError(f"unknown scenario '{scenario}'")
    base.update(overrides)
    return ScenarioConfig(**base)


# --------------------------------------------------------------------------- #
# 4) Per-subject generation.
# --------------------------------------------------------------------------- #
def _smooth_spatial(noise_ct: np.ndarray, coords: np.ndarray, sigma: float) -> np.ndarray:
    """Apply a spatial Gaussian smoother across channels (columns=time kept)."""
    if sigma <= 0:
        return noise_ct
    d2 = np.sum((coords[:, None, :] - coords[None, :, :]) ** 2, axis=-1)
    K = np.exp(-0.5 * d2 / sigma ** 2)
    K /= K.sum(axis=1, keepdims=True)
    return K @ noise_ct


def _individual_field(coords, times, rng, amp, rank):
    """A smooth low-rank per-subject spatiotemporal field = normal individual
    variability (NOT a quality problem)."""
    field = np.zeros((coords.shape[0], times.shape[0]))
    for _ in range(rank):
        fx = rng.uniform(-0.8, 0.8)
        fy = rng.uniform(-0.95, 0.95)
        spread = rng.uniform(0.35, 0.7)
        peak = rng.uniform(times[0], times[-1])
        width = rng.uniform(40, 120)
        sign = rng.choice([-1.0, 1.0])
        s = spatial_bump(coords, (fx, fy), spread)
        t = temporal_gaussian(times, peak, width)
        field += sign * np.outer(s, t)
    # normalize to unit peak, scale by amp
    peakabs = np.max(np.abs(field)) + 1e-9
    return amp * field / peakabs


def effect_map(cfg: ScenarioConfig, coords, times) -> np.ndarray:
    """The injected group-mean difference for beta 0: [nChan, nTime]."""
    if cfg.effect_amp == 0.0:
        return np.zeros((coords.shape[0], times.shape[0]))
    c1 = Component("eff1", cfg.effect_focus, cfg.effect_spread,
                   cfg.effect_peak_ms, cfg.effect_width_ms, cfg.effect_amp)
    c2 = Component("eff2", cfg.effect_focus2, cfg.effect_spread,
                   cfg.effect_peak_ms, cfg.effect_width_ms, cfg.effect_amp)
    return c1.render(coords, times) + c2.render(coords, times)


def simulate(cfg: ScenarioConfig):
    rng = np.random.default_rng(cfg.seed)
    names, coords = default_montage()
    nChan = coords.shape[0]
    times = np.linspace(cfg.t_start_ms, cfg.t_end_ms, cfg.nTime)
    adjacency = build_adjacency(coords)

    templates = base_template(coords, times, cfg.nBeta)   # [nBeta,nChan,nTime]
    eff = effect_map(cfg, coords, times)                  # [nChan,nTime], beta 0
    lib = _component_library()

    nSubject = cfg.n_per_group * cfg.n_groups
    groups = np.repeat(np.arange(cfg.n_groups), cfg.n_per_group).astype(np.int64)

    datain = np.zeros((cfg.nBeta, cfg.nTime, nChan, nSubject), dtype=np.float32)

    # ground-truth bookkeeping
    contaminated = np.zeros(nSubject, dtype=bool)
    contam_types = np.array(["none"] * nSubject, dtype=object)
    bad_channels = [[] for _ in range(nSubject)]
    sme = np.zeros(nSubject, dtype=np.float64)
    n_trials = np.zeros(nSubject, dtype=np.int64)
    latency = np.zeros(nSubject, dtype=np.float64)

    # decide which subjects get contaminated
    contam_pool = np.arange(nSubject)
    if cfg.contam_which_group >= 0:
        contam_pool = contam_pool[groups == cfg.contam_which_group]
    n_contam = int(round(cfg.contam_frac * len(contam_pool)))
    contam_ids = set(rng.choice(contam_pool, size=n_contam, replace=False).tolist()) \
        if (cfg.contam_type != "none" and n_contam > 0) else set()

    for s in range(nSubject):
        g = groups[s]
        # --- per-subject SIGNAL variability (must be preserved) ---
        lat = rng.normal(0.0, cfg.latency_sd_ms)
        # in latency_shift scenario, give group 1 a systematic latency offset
        if cfg.scenario == "latency_shift" and g == 1:
            lat += 30.0
        latency[s] = lat

        # rebuild templates with per-subject amplitude jitter + latency shift
        subj_beta = np.zeros((cfg.nBeta, nChan, cfg.nTime), dtype=np.float64)
        recipes = [["P1", "N170l", "N170r", "P3"], ["P1", "N2f", "P3"]]
        for b in range(cfg.nBeta):
            recipe = recipes[b] if b < len(recipes) else recipes[0]
            for cname in recipe:
                comp = lib[cname]
                amp_jit = comp.amp * (1.0 + rng.normal(0.0, cfg.amp_cv))
                jc = Component(comp.name, comp.focus, comp.spread,
                               comp.peak_ms, comp.width_ms, amp_jit)
                subj_beta[b] += jc.render(coords, times, latency_shift_ms=lat)
            # smooth individual field (normal individual differences)
            subj_beta[b] += _individual_field(coords, times, rng,
                                               cfg.indiv_field_amp, cfg.indiv_field_rank)

        # --- injected GROUP EFFECT (beta 0, groups >= 1) ---
        if g >= 1:
            subj_beta[0] += eff

        # --- measurement noise -> SME ---
        nt = int(max(10, rng.normal(cfg.n_trials_mean, cfg.n_trials_sd)))
        subj_sme = cfg.single_trial_sd / np.sqrt(nt)
        # group-level noise manipulation (whole group noisier)
        if g == 1:
            subj_sme *= cfg.group_noise_gain

        # --- CONTAMINATION (quality problems, should be down-weighted) ---
        ctype = "none"
        noise_gain = np.ones(nChan)
        if s in contam_ids:
            contaminated[s] = True
            ctype = cfg.contam_type
            if ctype == "variance_inflation":
                subj_sme *= cfg.var_inflation_gain
            elif ctype == "bad_channel":
                bad = rng.choice(nChan, size=cfg.bad_channel_k, replace=False)
                bad_channels[s] = [names[i] for i in bad]
                noise_gain[bad] = cfg.bad_channel_gain
                # also add a DC offset + drift to the bad channels
                for bch in bad:
                    drift = rng.normal(0, 3.0) * np.linspace(-1, 1, cfg.nTime)
                    subj_beta[:, bch, :] += drift[None, :]
            elif ctype == "time_burst":
                # transient artifact across a random channel subset
                k = rng.integers(4, 10)
                chs = rng.choice(nChan, size=k, replace=False)
                tprofile = temporal_gaussian(times, cfg.burst_peak_ms, cfg.burst_width_ms)
                for ch in chs:
                    burst = cfg.burst_amp * tprofile * rng.normal(1.0, 0.3)
                    subj_beta[:, ch, :] += burst[None, :]

        contam_types[s] = ctype
        sme[s] = subj_sme
        n_trials[s] = nt

        # add measurement noise (per beta), optionally spatially smoothed
        for b in range(cfg.nBeta):
            wn = rng.normal(0.0, 1.0, size=(nChan, cfg.nTime))
            wn = _smooth_spatial(wn, coords, cfg.noise_spatial_smooth)
            wn *= (subj_sme * noise_gain)[:, None]
            subj_beta[b] += wn

        # store as [nBeta, nTime, nChan]
        datain[:, :, :, s] = np.transpose(subj_beta, (0, 2, 1)).astype(np.float32)

    # effect mask (beta 0) from the injected map
    if eff.any():
        thr = cfg.effect_mask_frac * np.max(np.abs(eff))
        effect_mask = (np.abs(eff) >= thr)
    else:
        effect_mask = np.zeros((nChan, cfg.nTime), dtype=bool)

    out = dict(
        datain=datain,                      # [nBeta,nTime,nChan,nSubject]
        groups=groups,                      # [nSubject]
        adjacency=adjacency,                # [nChan,nChan]
        coords=coords,                      # [nChan,2]
        ch_names=np.array(names, dtype=object),
        times=times,                        # [nTime] ms
        effect_map=eff.astype(np.float32),  # [nChan,nTime]
        effect_mask=effect_mask,            # [nChan,nTime] bool (beta 0)
        contaminated=contaminated,          # [nSubject] bool (ground truth)
        contam_types=contam_types,          # [nSubject] str
        bad_channels=np.array(bad_channels, dtype=object),
        sme=sme,                            # [nSubject] measurement error
        n_trials=n_trials,                  # [nSubject]
        latency=latency,                    # [nSubject] ms shift
        config=cfg,
    )
    return out


# --------------------------------------------------------------------------- #
# 5) Saving.
# --------------------------------------------------------------------------- #
def save_outputs(out: dict, outdir: str):
    os.makedirs(outdir, exist_ok=True)
    cfg: ScenarioConfig = out["config"]

    # .npz for Python consumers
    npz_path = os.path.join(outdir, "sim_betas.npz")
    np.savez_compressed(
        npz_path,
        datain=out["datain"], groups=out["groups"], adjacency=out["adjacency"],
        coords=out["coords"], ch_names=out["ch_names"], times=out["times"],
        effect_map=out["effect_map"], effect_mask=out["effect_mask"],
        contaminated=out["contaminated"], contam_types=out["contam_types"].astype(str),
        sme=out["sme"], n_trials=out["n_trials"], latency=out["latency"],
    )

    # .mat for MATLAB / the LIMO side (optional but handy)
    mat_path = os.path.join(outdir, "sim_betas.mat")
    try:
        from scipy.io import savemat
        savemat(mat_path, {
            "datain": out["datain"], "groups": out["groups"] + 1,  # 1-based for MATLAB
            "adjacency": out["adjacency"], "coords": out["coords"],
            "ch_names": list(out["ch_names"]), "times": out["times"],
            "effect_mask": out["effect_mask"].astype(np.uint8),
            "contaminated": out["contaminated"].astype(np.uint8),
            "sme": out["sme"], "n_trials": out["n_trials"], "latency": out["latency"],
        }, do_compression=True)
    except Exception as e:  # scipy always present in limo-gae, but stay safe
        mat_path = f"(skipped: {e})"

    # human-readable manifest
    manifest = dict(
        scenario=cfg.scenario,
        shape_datain=list(out["datain"].shape),
        dims="[nBeta, nTime, nChan, nSubject]",
        n_groups=cfg.n_groups, n_per_group=cfg.n_per_group,
        nSubject=int(out["groups"].shape[0]),
        n_contaminated=int(out["contaminated"].sum()),
        contam_type=cfg.contam_type,
        effect_amp_uV=cfg.effect_amp,
        effect_peak_ms=cfg.effect_peak_ms,
        effect_cluster_size=int(out["effect_mask"].sum()),
        sme_mean=float(out["sme"].mean()), sme_min=float(out["sme"].min()),
        sme_max=float(out["sme"].max()),
        seed=cfg.seed,
        config=asdict(cfg),
    )
    with open(os.path.join(outdir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, default=str)

    return npz_path, mat_path, manifest


# --------------------------------------------------------------------------- #
# 6) CLI.
# --------------------------------------------------------------------------- #
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scenario", default="group_effect",
                   choices=["null", "group_effect", "bad_channel", "time_burst",
                            "variance_inflation", "latency_shift", "group_noise", "mixed"])
    p.add_argument("--n-per-group", type=int, default=20)
    p.add_argument("--n-groups", type=int, default=2)
    p.add_argument("--nBeta", type=int, default=2)
    p.add_argument("--effect-amp", type=float, default=None,
                   help="override injected effect amplitude (uV)")
    p.add_argument("--contam-frac", type=float, default=None)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default=None, help="output dir (default simulation/out/<scenario>)")
    args = p.parse_args()

    overrides = dict(n_per_group=args.n_per_group, n_groups=args.n_groups,
                     nBeta=args.nBeta, seed=args.seed)
    if args.effect_amp is not None:
        overrides["effect_amp"] = args.effect_amp
    if args.contam_frac is not None:
        overrides["contam_frac"] = args.contam_frac

    cfg = preset(args.scenario, **overrides)
    out = simulate(cfg)
    outdir = args.out or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                      "out", args.scenario)
    npz_path, mat_path, manifest = save_outputs(out, outdir)

    print(f"[simulate_group_betas] scenario = {cfg.scenario}")
    print(f"  datain shape [nBeta,nTime,nChan,nSubject] = {out['datain'].shape}")
    print(f"  groups: {cfg.n_groups} x {cfg.n_per_group} = {out['groups'].shape[0]} subjects")
    print(f"  contaminated subjects: {int(out['contaminated'].sum())} ({cfg.contam_type})")
    print(f"  injected effect: amp={cfg.effect_amp} uV @ {cfg.effect_peak_ms} ms, "
          f"cluster size={int(out['effect_mask'].sum())} (chan x time)")
    print(f"  SME range: [{out['sme'].min():.3f}, {out['sme'].max():.3f}] uV")
    print(f"  saved: {npz_path}")
    print(f"         {mat_path}")


if __name__ == "__main__":
    main()
