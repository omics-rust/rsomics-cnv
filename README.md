# rsomics-cnv

`rsomics-cnv` consolidates B-allele-frequency and log-R-ratio copy-number
workflows into one Rust product:

```text
rsomics-cnv call --sample TUMOR --control NORMAL --output calls input.bcf
rsomics-cnv polysomy --sample TUMOR --output polysomy input.bcf
```

Both commands accept VCF, BGZF-compressed VCF, BCF, and BGZF-compressed BCF.
They fail on malformed signals, unsorted records, ambiguous sample selection,
invalid model parameters, and an existing output directory. Report bundles are
committed as a complete directory and include a versioned JSON result alongside
the compatibility tables. `--json` emits the shared rsomics command envelope to
standard output.

The crate remains unpublished while region and target selection, external
allele-frequency input, plotting, and representative performance measurements
are completed. Unimplemented operations are not exposed by the CLI.

The compatibility reference is bcftools 1.24 `cnv`, `polysomy`, `HMM.c`, and
`peakfit.c`, retained under the upstream MIT license. Historical rsomics
implementations are used only as refactoring and fixture seeds.
