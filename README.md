# progressive-msa

A small, dependency-free Python library for progressive multiple sequence alignment (MSA), following the same overall strategy as Clustal Omega:

1. **Pairwise distances** via global Needleman-Wunsch alignment.
2. **Guide tree** construction via UPGMA clustering.
3. **Progressive alignment** of sequences/profiles by traversing the guide tree.

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

Pass `verbose=False` to `ProgressiveMSA(...)` to suppress progress logging.

## Development

```bash
pip install -e .
pip install pytest
pytest
```

## License

MIT — see [LICENSE](LICENSE).
