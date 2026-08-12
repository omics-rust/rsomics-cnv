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
the compatibility tables. The `call-result/v2` JSON is a compact artifact
manifest with chromosome and region summaries; per-site measurements and
posteriors remain in the compatibility tables instead of being duplicated in
JSON. `--json` emits the shared rsomics command envelope to standard output.
`call --allele-frequencies` accepts plain or gzip-compressed `CHROM`, `POS`,
`REF,ALT`, `AF` tables and restricts inference to listed sites.

Both operations accept inline or file-backed regions and targets. Regions use
TBI or CSI index jumps over BGZF VCF or BCF; targets stream over every accepted
input encoding and can be inverted with `^`. Position, record-span, and
variant-span overlap policies match bcftools 1.24. BED, VCF, and generic
tabular files use the shared `rsomics-intervals` coordinate contract.

`call --optimize FRACTION` estimates chromosome-specific aberrant-cell
fractions and BAF deviations for query and control samples, records them in
compatibility and JSON reports, and falls back to the declared starting model
when the iterative fit does not converge. `-a`, `-d`, and `-k` accept either a
query value or `QUERY,CONTROL` values for matched samples.

`call -p QUALITY` adds LRR, BAF, and copy-number SVGs for chromosomes reaching
the requested region quality. `polysomy --plots` adds fitted BAF-distribution
SVGs and a chromosome copy-number overview. The plots are self-contained and
do not require Python or a plotting runtime.

The CLI infers, writes, and releases one chromosome at a time. On the tracked
300,000-site macOS arm64 release gate, `call` used 32.1% less peak memory than
bcftools 1.24 and `polysomy` was 4.4 times faster while using 49.7% less peak
memory. These are local equivalent-workflow measurements, not cross-platform
claims; commands, raw measurements, fixture hashes, and compatibility checks
are retained in [`benchmarks`](benchmarks/README.md).

Unimplemented behavior is not exposed by the CLI. Every release requires a
fresh package review and exact-head four-platform CI.

The compatibility reference is bcftools 1.24 `cnv`, `polysomy`, `HMM.c`, and
`peakfit.c`, retained under the upstream MIT license. Historical rsomics
implementations are used only as refactoring and fixture seeds.
