# Tests

Run `make check` for the active joint-inference regression suite in `joint/`.
It covers exact toy rewards and certificates, paired molecule semantics,
overflow guards, bounded descriptor use, serial/parallel equivalence, fresh
SR/Hi-C/mixed inputs on DAG and cyclic graphs, graph/index round trips,
multi-component binning and unsupported-option rejection.

The raw integration tests generate their own sequences and paired FASTQs.
They verify emitted haplotypes and independently replay the objective.
They need no downloaded benchmark data and no archived inferred solution.

The full raw-read-to-VCF path, with the embedded caller and an empty `PATH`:

```bash
python3 test/joint/raw_toy_test.py --executable ./DipGenie --libraries mixed --recombinant --dipcall
```

The MHC files in this directory are the inputs of the test run in the main
[README](../README.md): a four-haplotype MHC panel (`MHC_4.vcf.gz`, `MHC_4.gfa.gz`),
its reference (`MHC-CHM13.0.fa.gz`) and HG002 paired short reads at 2x
(`HG002.mhc.R1.fq.gz`, `HG002.mhc.R2.fq.gz`).
