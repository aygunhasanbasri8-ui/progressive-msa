"""Command line entry point: `msa align input.fasta -o out.aln --format clustal`."""

from __future__ import annotations

import argparse
import logging
import sys

from .evaluate import sum_of_pairs, total_column
from .io import FastaError, format_clustal, format_fasta, parse_fasta, read_fasta
from .msa_core import ProgressiveMSA
from .scoring import Scoring

FORMATTERS = {"fasta": format_fasta, "clustal": format_clustal}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="msa", description="Progressive multiple sequence alignment.")
    sub = parser.add_subparsers(dest="command", required=True)

    align = sub.add_parser("align", help="align the sequences of a FASTA file")
    align.add_argument("input", help="FASTA file, or '-' for stdin")
    align.add_argument("-o", "--output", help="output file (default: stdout)")
    align.add_argument("-f", "--format", choices=sorted(FORMATTERS), default="fasta")
    align.add_argument(
        "--tree-order",
        action="store_true",
        help="list sequences in guide-tree order instead of input order",
    )
    align.add_argument(
        "-s",
        "--scoring",
        choices=("unit", "blosum62", "dna"),
        default="unit",
        help="unit: edit distance (default); blosum62: proteins; dna: transitions cheaper than transversions",
    )
    align.add_argument("--gap-open", type=float, default=10.0, help="gap opening penalty (blosum62/dna)")
    align.add_argument("--gap-extend", type=float, default=1.0, help="gap extension penalty (blosum62/dna)")
    align.add_argument("-t", "--tree", choices=("upgma", "nj"), default="upgma", help="guide tree method")
    align.add_argument("-v", "--verbose", action="store_true", help="log the pipeline steps to stderr")

    score = sub.add_parser("score", help="compare an alignment with a reference (SP and TC scores)")
    score.add_argument("test", help="aligned FASTA to evaluate")
    score.add_argument("reference", help="reference alignment in aligned FASTA")
    return parser


def _align(args: argparse.Namespace) -> int:
    try:
        if args.input == "-":
            sequences = parse_fasta(sys.stdin)
        else:
            sequences = read_fasta(args.input)
    except (OSError, FastaError) as exc:
        print(f"msa: {exc}", file=sys.stderr)
        return 2

    scoring = None
    if args.scoring == "blosum62":
        scoring = Scoring.blosum62(args.gap_open, args.gap_extend)
    elif args.scoring == "dna":
        scoring = Scoring.dna(gap_open=args.gap_open, gap_extend=args.gap_extend)
    aligned = ProgressiveMSA(sequences, verbose=False, scoring=scoring, tree=args.tree).run()
    if not args.tree_order:
        aligned = {name: aligned[name] for name in sequences}

    text = FORMATTERS[args.format](aligned)
    if args.output:
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(text)
    else:
        sys.stdout.write(text)
    return 0


def _score(args: argparse.Namespace) -> int:
    try:
        test = read_fasta(args.test, keep_gaps=True)
        reference = read_fasta(args.reference, keep_gaps=True)
        sp, tc = sum_of_pairs(test, reference), total_column(test, reference)
    except (OSError, ValueError) as exc:
        print(f"msa: {exc}", file=sys.stderr)
        return 2
    print(f"SP {sp:.3f}  TC {tc:.3f}")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "score":
        return _score(args)
    if args.verbose:
        logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stderr)
    return _align(args)


if __name__ == "__main__":
    raise SystemExit(main())
