"""Alignment quality against a reference (BAliBASE-style SP and TC scores)."""

from __future__ import annotations

from collections.abc import Mapping


def _columns(alignment: Mapping[str, str]) -> list[frozenset[tuple[str, int]]]:
    """Each column as the set of (sequence, residue index) it aligns; gaps left out."""
    positions = dict.fromkeys(alignment, 0)
    length = len(next(iter(alignment.values()))) if alignment else 0
    columns = []
    for k in range(length):
        column = []
        for name, row in alignment.items():
            if row[k] != "-":
                column.append((name, positions[name]))
                positions[name] += 1
        columns.append(frozenset(column))
    return columns


def _pairs(columns) -> set[tuple[tuple[str, int], tuple[str, int]]]:
    pairs = set()
    for column in columns:
        items = sorted(column)
        for x, a in enumerate(items):
            for b in items[x + 1 :]:
                pairs.add((a, b))
    return pairs


def _check(test: Mapping[str, str], reference: Mapping[str, str]) -> None:
    if set(test) != set(reference):
        raise ValueError("test and reference must contain the same sequence names")
    for name in reference:
        if test[name].replace("-", "") != reference[name].replace("-", ""):
            raise ValueError(f"sequence {name!r} differs between test and reference")


def sum_of_pairs(test: Mapping[str, str], reference: Mapping[str, str]) -> float:
    """Share of residue pairs aligned in the reference that the test also aligns (0-1)."""
    _check(test, reference)
    ref = _pairs(_columns(reference))
    if not ref:
        return 1.0
    return len(ref & _pairs(_columns(test))) / len(ref)


def total_column(test: Mapping[str, str], reference: Mapping[str, str]) -> float:
    """Share of reference columns (with at least two residues) reproduced exactly (0-1)."""
    _check(test, reference)
    ref = [c for c in _columns(reference) if len(c) >= 2]
    if not ref:
        return 1.0
    found = set(_columns(test))
    return sum(c in found for c in ref) / len(ref)
