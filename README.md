# progressive-msa

English | [Türkçe](README.tr.md)

[![CI](https://github.com/aygunhasanbasri8-ui/progressive-msa/actions/workflows/ci.yml/badge.svg)](https://github.com/aygunhasanbasri8-ui/progressive-msa/actions/workflows/ci.yml)

A small, dependency-free Python library for progressive multiple sequence alignment (MSA), following the same overall strategy as Clustal Omega:

1. **Pairwise distances** via global Needleman-Wunsch alignment.
2. **Guide tree** construction via UPGMA clustering.
3. **Progressive alignment** of sequences/profiles by traversing the guide tree.

It is a compact, readable implementation of the classic pipeline, meant for learning and for
small inputs. See [Limitations](#limitations) before using it on real biological data.

## Features

- Pure Python, no runtime dependencies, Python 3.9 or newer.
- Deterministic results: ties in the alignment, in the consensus and in the UPGMA pair selection are always broken
  the same way, independent of `PYTHONHASHSEED`.
- Every step is also callable on its own (`calculate_distance_matrix`, `build_guide_tree`, `progressive_alignment`).
- Optional progress logging (printed with `verbose=True`, always sent to the `progressive_msa` logger).
- FASTA input, FASTA or Clustal output, and an `msa` command.

## Installation

```bash
pip install -e .
```

## Usage

```python
from progressive_msa import ProgressiveMSA

sequences = {
    "Seq1": "ACTG",
    "Seq2": "ACTC",
    "Seq3": "ACCGT",
}

msa = ProgressiveMSA(sequences)
aligned = msa.run()

for name, seq in aligned.items():
    print(name, seq)
```

Output (the progress log is printed first unless `verbose=False`):

```
Seq3 ACCGT
Seq1 ACTG-
Seq2 ACTC-
```

Pass `verbose=False` to `ProgressiveMSA(...)` to suppress progress logging.

The returned dictionary lists the sequences in the order they are visited in the guide tree,
not in input order.

### Command line

```bash
msa align input.fasta                       # aligned FASTA to stdout
msa align input.fasta -o out.aln -f clustal  # Clustal-style blocks with a conservation line
cat input.fasta | msa align -               # read from stdin
```

Sequences are written in input order; `--tree-order` keeps the guide-tree order instead.
Gaps in the input are ignored, the sequences are aligned from scratch. `-v` logs the pipeline
steps to stderr.

From Python: `parse_fasta`, `read_fasta`, `format_fasta` and `format_clustal`.

### Scoring

By default the library keeps its original model: unit costs and an edit distance. For real
sequences pass a scoring scheme; the pairwise step then uses Gotoh's affine-gap algorithm
(gap opening and extension penalised separately) and the guide tree is built from
`1 - identity` of those alignments:

```python
from progressive_msa import ProgressiveMSA, Scoring

ProgressiveMSA(proteins, scoring=Scoring.blosum62())  # gap open 10, extend 1
ProgressiveMSA(proteins, scoring=Scoring.blosum62(gap_open=11))
ProgressiveMSA(reads, scoring=Scoring.dna())  # match 5, transition -1, transversion -4
```

On the command line: `msa align in.fasta --scoring blosum62 --gap-open 10 --gap-extend 1`.
BLOSUM62 is the table NCBI distributes. Affine gaps keep an insertion or deletion in one
piece where unit costs would scatter it.

### API

| Method | Description |
|---|---|
| `ProgressiveMSA(sequences, verbose=True)` | `sequences` maps a name to a sequence string |
| `run()` | Runs all three stages and returns `{name: aligned_sequence}` |
| `calculate_distance_matrix()` | Pairwise edit-distance matrix (list of lists, symmetric, zero diagonal) |
| `needleman_wunsch_distance(seq1, seq2)` | Global-alignment edit distance of two sequences |
| `build_guide_tree()` | UPGMA guide tree as nested tuples, e.g. `('Seq3', ('Seq1', 'Seq2'))`; call after `calculate_distance_matrix()` |
| `progressive_alignment()` | Aligns along the guide tree; call after `build_guide_tree()` |

## How it works

1. **Distance matrix.** Every pair of sequences is aligned with Needleman-Wunsch using unit
   costs (mismatch 1, gap 1); the resulting edit distance fills a symmetric matrix.
2. **Guide tree.** UPGMA repeatedly merges the closest pair of clusters and averages their
   distances to the other clusters.
3. **Progressive alignment.** The tree is traversed bottom-up. At each internal node the two
   child profiles are aligned through their consensus strings (the most frequent residue
   per column), and the resulting gap pattern is applied to every sequence in both profiles.

## Limitations

- The default scoring is simple (mismatch and gap both cost 1); use `Scoring.blosum62()` or `Scoring.dna()` for substitution scores and affine gaps. PAM matrices are not included.
- With the default unit costs profiles are aligned through their consensus sequence; with a `Scoring` they are aligned column to column.
- The distance matrix takes O(n²·L²) time for n sequences of length L, so it is not meant for large inputs.
- UPGMA (default) assumes a constant rate of evolution; `tree="nj"` uses neighbor joining instead.

## Guide tree and profiles

`ProgressiveMSA(seqs, scoring=Scoring.blosum62(), tree="nj")` builds the guide tree with
neighbor joining (Saitou-Nei), which does not assume a constant rate of evolution. With a
scoring scheme, profiles are aligned column to column: two columns score the average
substitution score over every residue pair across them, so minority residues still count
(the default unit-cost mode keeps aligning consensus strings).

## Measuring quality

`msa score test.fa reference.fa` prints the sum-of-pairs (SP: share of reference residue
pairs the test also aligns) and total-column (TC: share of reference columns reproduced
exactly) scores, the measures used by BAliBASE. From Python: `sum_of_pairs` and
`total_column`.

`benchmarks/simulate.py` generates protein families with a known true alignment
(substitutions plus short insertions/deletions along random lineages) and scores each
configuration against it. With 12 families of 6 sequences, length 120, 25% substitutions
and 2% indels per residue:

| Configuration | SP | TC |
|---|---|---|
| unit costs / UPGMA (default) | 0.935 | 0.813 |
| BLOSUM62 + affine gaps, column profiles / UPGMA | 0.948 | 0.850 |
| BLOSUM62 + affine gaps, column profiles / NJ | 0.942 | 0.842 |

Substitutions in the simulation are uniform rather than BLOSUM-shaped, so this
understates what a substitution matrix gains on real proteins; at 45% substitutions the
gap is similar (SP 0.751 vs 0.768). A comparison on BAliBASE against Clustal Omega and
MUSCLE has not been done yet.

## Development

```bash
pip install -e ".[dev]"
ruff check . && ruff format --check .
pytest
```

Publishing a GitHub release `vX.Y.Z` builds the package and uploads it to PyPI through
Trusted Publishing (`.github/workflows/publish.yml`; the project must be registered on PyPI
with this workflow first).

GitHub Actions runs ruff and the tests on Python 3.9, 3.11 and 3.12 for every push and pull request
(`.github/workflows/ci.yml`).

## License

MIT. See [LICENSE](LICENSE).
