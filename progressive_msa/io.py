"""Reading and writing alignments: FASTA in, FASTA or Clustal out."""

from __future__ import annotations

from collections.abc import Iterable, Mapping

CLUSTAL_BLOCK = 60
FASTA_LINE = 60


class FastaError(ValueError):
    """Raised when the input is not usable FASTA."""


def parse_fasta(lines: Iterable[str]) -> dict[str, str]:
    """Parses FASTA records into {name: sequence}, keeping input order.

    The name is the first word of the header line, so "seq1 Homo sapiens" becomes "seq1".
    """
    sequences: dict[str, str] = {}
    name = None
    chunks: list[str] = []

    def flush() -> None:
        if name is None:
            return
        sequence = "".join(chunks).upper()
        if not sequence:
            raise FastaError(f"record '{name}' has no sequence")
        sequences[name] = sequence

    for line_no, raw in enumerate(lines, start=1):
        line = raw.strip()
        if not line or line.startswith(";"):
            continue
        if line.startswith(">"):
            flush()
            header = line[1:].split()
            if not header:
                raise FastaError(f"line {line_no}: empty header")
            name = header[0]
            if name in sequences:
                raise FastaError(f"line {line_no}: duplicate name '{name}'")
            chunks = []
        elif name is None:
            raise FastaError(f"line {line_no}: sequence data before the first '>' header")
        else:
            # Gaps from an already aligned input are dropped; the library aligns from scratch.
            chunks.append("".join(line.split()).replace("-", "").replace(".", ""))

    flush()
    if not sequences:
        raise FastaError("no FASTA records found")
    return sequences


def read_fasta(path: str) -> dict[str, str]:
    with open(path, encoding="utf-8") as handle:
        return parse_fasta(handle)


def format_fasta(alignment: Mapping[str, str], width: int = FASTA_LINE) -> str:
    out = []
    for name, seq in alignment.items():
        out.append(f">{name}")
        out.extend(seq[i : i + width] for i in range(0, len(seq), width))
    return "\n".join(out) + "\n"


def _conservation(column: tuple[str, ...]) -> str:
    first = column[0]
    return "*" if first != "-" and all(c == first for c in column) else " "


def format_clustal(alignment: Mapping[str, str], block: int = CLUSTAL_BLOCK) -> str:
    """Clustal-style output: blocks of `block` columns and a '*' line for identical columns."""
    if not alignment:
        return "CLUSTAL multiple sequence alignment (progressive-msa)\n"

    names = list(alignment)
    pad = max(len(n) for n in names) + 4
    length = len(next(iter(alignment.values())))
    out = ["CLUSTAL multiple sequence alignment (progressive-msa)", "", ""]

    for start in range(0, length, block):
        rows = {n: alignment[n][start : start + block] for n in names}
        out.extend(f"{n:<{pad}}{rows[n]}" for n in names)
        marks = "".join(_conservation(col) for col in zip(*rows.values()))
        out.append(" " * pad + marks)
        out.append("")

    return "\n".join(out) + "\n"
