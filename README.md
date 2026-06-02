# rsomics-phydiv

Generalized phylogenetic alpha-diversity from a feature count table and a rooted
Newick tree. Computes the McCoy & Matsen balance-weighted PD family (BWPD): the
rooted/unrooted choice and the abundance-weighting exponent θ are both exposed,
so a single binary covers Faith's PD, unrooted PD, and the weighted variants.

Faith's PD (rooted, unweighted) has its own crate, [`rsomics-faith-pd`]; this is
the superset operation.

## Usage

```sh
rsomics-phydiv counts.tsv --tree tree.nwk                       # rooted Faith's PD (auto)
rsomics-phydiv counts.tsv --tree tree.nwk --unrooted            # unrooted PD
rsomics-phydiv counts.tsv --tree tree.nwk --weight 1            # rooted, fully weighted
rsomics-phydiv counts.tsv --tree tree.nwk --unrooted --weight 0.25   # BWPD_0.25
```

The table is feature-by-sample TSV (first column feature IDs, header row sample
names). Output is `sample<TAB>phydiv`, one row per sample.

`--rooted` / `--unrooted` override the default, which follows scikit-bio: rooted
iff the tree's root is bifurcating. `--weight` is 0 (unweighted), 1 (fully
weighted), or a θ in (0, 1) for partial weighting.

## Origin

This crate is an independent Rust reimplementation of `scikit-bio`'s
`skbio.diversity.alpha.phydiv` based on:

- The published methods: Faith 1992 (DOI 10.1016/0006-3207(92)91201-3),
  McCoy & Matsen 2013 balance-weighted PD (DOI 10.7717/peerj.157), and the
  array-form node-abundance formulation of Hamady, Lozupone & Knight 2010
  (DOI 10.1038/ismej.2009.97).
- The scikit-bio source (BSD-3-Clause), which was read and cited for the exact
  rooted-vs-unrooted branch selection, the `2·min(p, 1−p)` balance fold, the
  θ-power semantics (skipped at θ=0 and θ=1), and the missing-length-as-zero
  convention.

Compatibility is value-exact against scikit-bio 0.7.2 within 1e-9 (deterministic
branch-length summation weighted by closed-form node abundances).

License: MIT OR Apache-2.0.
Upstream credit: [scikit-bio](https://github.com/scikit-bio/scikit-bio)
(BSD-3-Clause).

[`rsomics-faith-pd`]: https://github.com/omics-rust/rsomics-faith-pd
