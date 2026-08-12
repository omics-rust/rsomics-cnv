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
standard output. `call --allele-frequencies` accepts plain or gzip-compressed
`CHROM`, `POS`, `REF,ALT`, `AF` tables and restricts inference to listed sites.

Both operations accept inline or file-backed regions and targets. Regions use
TBI or CSI index jumps over BGZF VCF or BCF; targets stream over every accepted
input encoding and can be inverted with `^`. Position, record-span, and
variant-span overlap policies match bcftools 1.24. BED, VCF, and generic
tabular files use the shared `rsomics-intervals` coordinate contract.

The crate remains unpublished while call optimization, plotting, and
representative equivalent-workflow timing and peak-memory measurements are
completed. Unimplemented behavior is not exposed by the CLI.

The compatibility reference is bcftools 1.24 `cnv`, `polysomy`, `HMM.c`, and
`peakfit.c`, retained under the upstream MIT license. Historical rsomics
implementations are used only as refactoring and fixture seeds.
