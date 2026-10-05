# limo_tfce tests

Self-contained validation of the exact eTFCE engine in `limo_tfce.m`
(no Image Processing Toolbox / SPM required). Run from MATLAB with the repo
on the path:

```matlab
test_limo_etfce        % exact TFCE (minnbchan = 0) vs an independent
                       % brute-force exact integrator (graph/conncomp),
                       % for type 1/2/3, signed and positive maps, and stacks
test_minnbchan_etfce   % minnbchan pruning (k-core "entry height" h_eff) vs an
                       % independent exact reference that prunes with the direct
                       % iterative FieldTrip rule; checks default == minnbchan=2
```

Both compare the union-find forest result against a completely independent
implementation and pass at ~1e-15 (machine precision). They do not exercise
`method = 'classic'` (that path is the original, unchanged code and needs
`limo_findcluster`).
