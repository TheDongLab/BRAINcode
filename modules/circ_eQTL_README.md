# Optional circRNA abundance mode

From the QTL scripts directory:

```bash
sbatch prep_eQTL.sh Cerebellum circ
# After preparation completes successfully:
sbatch run_eQTL.sh Cerebellum circ
# Optional interaction analysis, after the standard run:
sbatch run_eQTL.sh Cerebellum circ interaction
```

Commands without `circ` retain the ordinary gene-expression workflow.
The runner also accepts `circ` after the legacy run-type/subdirectory arguments.
Circ-mode inputs and results live under `<tissue>/circ_eQTL/`.

Preparation reads the existing `QTL/circ_rpm_matrix.txt`, applies tissue remapping
(including the observed `Lumbosacaral` spelling), RIN >=3, PMI <=40, sex mismatch
exclusions, and highest-RIN sample selection per donor. It intersects the selected
samples with genotype and complete covariate data. The final donor list is written
to `selected_samples.tsv`.

There is **no minimum BSJ-read or prevalence threshold**. Only all-zero/constant
RPM rows are removed. RPM is rank inverse-normal transformed per circRNA across
the final tissue cohort, using average ranks for ties and `(rank-0.5)/N`.
`expression_<tissue>.txt` contains the transformed phenotype; `circ_rpm_<tissue>.txt`
preserves the same rows/donors in original RPM units. Plots use the transformed
phenotype and label it accordingly. Raw BSJ counts remain in the master count matrix.

Coordinates are parsed from circRNA IDs. The original BED start is converted to
one-based `start+1`; end is unchanged. The existing 1-Mb cis window applies around
this circRNA interval. IDs remain unchanged. `circ_annotation.tsv` contains host
genes/transcripts collected from the original CIRCexplorer2 files listed in
`circ_rpm_sample_metadata.tsv` for selected samples. Those files must remain accessible.

Circ mode uses the existing SNP/covariate/model settings. SNP MAF and genotype-group
filters are evaluated after the complete-covariate donor intersection. Matrices
must have exactly matching donor IDs/order. Circ FDR outputs contain only actual
FDR-passing associations; the legacy nominal-p-value fallback is not used in circ mode.
The regional plot retains the existing predefined chr17 window and skips it cleanly
if no circ associations occur there.

Tests: `python3 tests/test_circ_eqtl.py`, `python3 tests/test_eqtl_modes.py`, and
`Rscript tests/smoke_circ_eqtl.R modules` (requires the pipeline's Python/R packages).
