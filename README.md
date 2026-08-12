# rsomics-cnv

`rsomics-cnv` consolidates B-allele-frequency and log-R-ratio copy-number
workflows into one Rust product. The target command family is:

```text
rsomics-cnv call
rsomics-cnv polysomy
```

The crate is under active reconstruction and is not published. The binary
will remain absent until both commands have complete typed input, checked
models, transactional reports, bcftools 1.24 compatibility evidence, and
representative performance measurements.

The compatibility reference is bcftools 1.24 `cnv`, `polysomy`, `HMM.c`, and
`peakfit.c`, retained under the upstream MIT license. Historical rsomics
implementations are used only as refactoring and fixture seeds.
