# Directional Coherence in Heterogeneous Cosmological Anomalies

Verified catalogue, consolidated analysis code, and pre-registration for the
cosmological companion paper to *Nothing, Described: A Structural Ontology of
Everything* (Biru et al. 2026a, 2026b).

## Contents

    data/verified_independent.csv   The citable dataset: 6 verified independent
                                    anomalies + 1 contested (GRB H0 dipole),
                                    every coordinate checked against its
                                    primary measurement paper (August 2026).
    data/catalogue_annotated.csv    Wider catalogue of entries handled by the
                                    programme, each tagged with status:
                                    core / contested / excluded_* /
                                    collapsed_* / removed_unsourced. Unverified
                                    coordinates are marked as approximate.
    code/analysis.py                One script reproducing every statistic in
                                    the companion paper.
    preregistration/                The OSF pre-registration document.
    results/                        Generated: summary + figures.

## Reproducing the paper

    pip install -r requirements.txt
    python code/analysis.py          # full Monte Carlo (~minutes)
    python code/analysis.py --fast   # reduced Monte Carlo, quick check

The script prints and writes: directional clustering (Monte Carlo),
K_simple clustering, the fanning correlation with exact permutation
p-values, the full-sky reference scan with look-elsewhere correction,
functional-form fits with AIC, the trunk extrapolation with full
covariance propagation, and the four-point additive-tree condition with
matched-redshift random baselines. RNG is seeded (seed = 7); Monte Carlo
p-values match the paper within sampling noise.

## Provenance and honesty notes

The catalogue underwent a verification audit in August 2026 (companion
paper, Section 3.4): two coordinates corrected, one entry removed as
unsourced, one entry reclassified as contested. Earlier statistics
computed on pre-verification catalogues are superseded and are retained
only in the papers' audit sections, deliberately, as the record of the
correction. The removed entry appears in catalogue_annotated.csv with
status removed_unsourced so that the audit trail is public.

## Licence

CC BY 4.0. See LICENCE and CITATION.cff.
