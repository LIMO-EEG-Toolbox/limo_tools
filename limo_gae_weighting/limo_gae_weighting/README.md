# Bounded group-blind quality weights for LIMO's second level

```
gae_weighting.py            the method + the estimators
demo_homogeneous_split.py   one homogeneous cohort -> random split -> inject -> compare
simulate_group_betas.py     synthetic beta-map generator (32-channel montage)
```

Method: masked-channel graph autoencoder on montage adjacency (k = 6), out-of-fold
reconstruction error per (channel, participant), robust z within channel (median/MAD),
one-sided Tukey bisquare with tau = 3, c = 4.685, mean 1 within channel. Output `W[nChan, nSub]`
goes into LIMO's weighted second level. Group labels are never used.

```bash
pip install numpy scipy torch torch_geometric

# weights for your own data (npz with `datain[nBeta,nTime,nChan,nSub]` + `adjacency` or `coords`)
python gae_weighting.py my_betas.npz --beta 0 --out weights.npz

# validation protocol
python demo_homogeneous_split.py                                          # inert?
python demo_homogeneous_split.py --effect-amp 2.5                         # effect retained?
python demo_homogeneous_split.py --contam-frac 0.3 --contam-kind offset   # fabrication?
```

The demo prints JSON: flag rate per arm (Mann-Whitney), hit rate on corrupted cells, effect
retained, significant cells, permutation FWER p, for mean / Yuen 20 % trim / weighted.

Open question for the student: the weights are label-free but data-value-dependent, so does the
weighted Welch t (Kish effective sample size) keep its nominal null? Axes: n per arm,
contamination prevalence, unequal variances, heavy tails, effect shape (including a latency
shift, which must not be down-weighted: `simulate_group_betas.py --scenario latency_shift`),
adjacency k (`--knn`), tau.
