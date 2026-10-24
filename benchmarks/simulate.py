"""Simulated families with a known true alignment, scored with SP/TC.

Each family starts from a random protein; new members are copied from an existing
member and mutated (substitutions, short insertions and deletions) while the true
alignment is tracked column by column. The aligner never sees it.

    python benchmarks/simulate.py --families 20 --members 6 --length 120
"""

from __future__ import annotations

import argparse
import random
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from progressive_msa import ProgressiveMSA, Scoring, sum_of_pairs, total_column  # noqa: E402

AMINO = "ARNDCQEGHILKMFPSTWYV"


def simulate(rnd: random.Random, members: int, length: int, sub_rate: float, indel_rate: float):
    # rows[name] is a list over the shared columns: residue or "-"
    rows = {"s0": [rnd.choice(AMINO) for _ in range(length)]}
    for k in range(1, members):
        parent = rows[rnd.choice(list(rows))]
        child = [c if c == "-" or rnd.random() > sub_rate else rnd.choice(AMINO) for c in parent]
        rows[f"s{k}"] = child
        col = 0
        while col < len(child):
            if child[col] != "-" and rnd.random() < indel_rate:
                size = rnd.randint(1, 4)
                if rnd.random() < 0.5:
                    for c in range(col, min(col + size, len(child))):
                        child[c] = "-"
                else:
                    for row in rows.values():
                        for _ in range(size):
                            row.insert(col + 1, rnd.choice(AMINO) if row is child else "-")
                    col += size
            col += 1
    keep = [c for c in range(len(rows["s0"])) if any(row[c] != "-" for row in rows.values())]
    reference = {name: "".join(row[c] for c in keep) for name, row in rows.items()}
    sequences = {name: aligned.replace("-", "") for name, aligned in reference.items()}
    return sequences, reference


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--families", type=int, default=20)
    p.add_argument("--members", type=int, default=6)
    p.add_argument("--length", type=int, default=120)
    p.add_argument("--sub-rate", type=float, default=0.25)
    p.add_argument("--indel-rate", type=float, default=0.02)
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()

    configs = {
        "unit / UPGMA (default)": dict(),
        "BLOSUM62 + affine / UPGMA": dict(scoring=Scoring.blosum62()),
        "BLOSUM62 + affine / NJ": dict(scoring=Scoring.blosum62(), tree="nj"),
    }
    rnd = random.Random(args.seed)
    families = [
        simulate(rnd, args.members, args.length, args.sub_rate, args.indel_rate) for _ in range(args.families)
    ]

    print(
        f"{args.families} families x {args.members} members, length {args.length}, "
        f"substitution {args.sub_rate}, indel {args.indel_rate} per residue"
    )
    print(f"{'configuration':30} {'SP':>6} {'TC':>6} {'s/family':>9}")
    for label, options in configs.items():
        sp, tc = [], []
        start = time.perf_counter()
        for sequences, reference in families:
            aligned = ProgressiveMSA(sequences, verbose=False, **options).run()
            sp.append(sum_of_pairs(aligned, reference))
            tc.append(total_column(aligned, reference))
        per = (time.perf_counter() - start) / len(families)
        print(f"{label:30} {statistics.mean(sp):6.3f} {statistics.mean(tc):6.3f} {per:9.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
