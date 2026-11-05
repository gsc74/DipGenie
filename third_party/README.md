# Third-party code

Initialize pinned dependencies with `git submodule update --init --recursive`.

- `minimap2` v2.31: assembly alignment for the embedded dipcall stages. Its
  static library is linked. Private kseq
  symbols are prefixed to avoid clashes.
- `htsbox/htslib` (HTSbox r352): indexed VCF access and BAM I/O. The pinned HTSbox pileup,
  FASTA indexing and RAZF sources are linked for dipcall assembly calling;
  `src/EmbeddedPileup.c` redirects only the library's stream output and
  serializes its getopt entry point, without process-wide stdout redirection.
- `src/AssemblyCalling.cpp`: dipcall 0.3 samflt/vcfpair policies adapted from
  `6bd5d7724699491f215aeb5fb628490ebf2cc3ae`; MIT notice is in
  `DIPCALL-LICENSE.txt`. Callable regions use the pinned minimap2 paftools policy.
- `unordered_dense` v5.3.3: header-only hash maps (MIT).
- `src/PanGenieVcf.*`: phased panel import adapted from PanGenie, with its MIT notice.
- `src/gfa-base.cpp` and supporting headers: graph construction utilities
  derived from minigraph; notices are retained.
- `src/NativeLocalizer.hpp`: narrow exact seed-extension and interval-voting
  adaptation from hifiasm `ec9a8b222d149d25b7355e83765698640d59b189` (MIT).
  Physical graph coordinates and selected-anchor APIs are DipGenie adapters;
  no full assembler dependency is imported. See `HIFIASM-LICENSE.txt`.

The inference engine contains no alternative caller or long-read inference
path, and read localization invokes no external mapper. `VariantOutput.cpp` handles
the internal caller's literal-allele GT/AD VCF contract, final phase-set fields,
component joining, BGZF and CSI. It is differential-tested against bcftools,
which is a test oracle, not a runtime dependency. PAV and Snakemake are
not dependencies.

See the root LICENSE and the retained notices in this directory and each
submodule for licensing terms.
