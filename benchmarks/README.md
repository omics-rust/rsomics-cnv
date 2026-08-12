# Release performance gate

The release workload contains 12 chromosomes and 300,000 sites per workflow.
`call` uses paired `QUERY` and `CONTROL` samples with CN1, CN2, and CN3
segments and chromosome-level aberrant-fraction optimization. `polysomy`
contains CN2, CN3, and fractional CN3.5 decisions.

The fixture generator is `generate_release_fixtures.awk`. The tracked raw
measurements are from 12 paired rounds with alternating execution order:

```text
benchmarks/run_release_gate.sh FIXTURE_DIRECTORY RESULT_TSV 12
```

Each round compares `call` BAF/LRR and posterior tables byte for byte, then
checks cell-fraction estimates and region coordinates and states. For
`polysomy`, it compares every normalized distribution bin and copy-number
decision, with a `1e-5` tolerance for fitted deviations. A failed comparison
aborts the measurement.

## 2026-08-13 macOS arm64

Apple M2, macOS 26.6.1, bcftools 1.24, release profile. Times are seconds and
RSS is the median maximum resident set size reported by `/usr/bin/time -lp`.

| Workflow | Tool | Wall median | CPU median | RSS median |
|---|---:|---:|---:|---:|
| call | bcftools | 1.375 | 1.200 | 37,978,112 B |
| call | rsomics-cnv | 1.420 | 1.170 | 25,772,032 B |
| polysomy | bcftools | 5.160 | 5.025 | 7,536,640 B |
| polysomy | rsomics-cnv | 1.170 | 1.090 | 3,792,896 B |

`call` reduces peak memory by 32.1% and CPU time by 2.5%; its median wall time
is 3.3% higher. `polysomy` reduces wall time by 77.3% and peak memory by 49.7%.
The release decision rests on strict peak-memory improvement for `call` and on
both throughput and memory improvements for `polysomy`.

Raw samples and complete provenance are in
`results/2026-08-13-macos-aarch64.tsv` and its `.meta` companion. The fixture
SHA-256 values are:

```text
fb5b1fa37ea19d60fafe12320a5fa8edfb18c65e732f46b0a080d3f9b8a4eb76  call.vcf
b15942724f7e47c8648830125400b5bf4c2eb1f446a09f863e4af71a358b8c1d  polysomy.vcf
```
