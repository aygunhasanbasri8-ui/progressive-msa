"""Scoring schemes: substitution matrices and affine gap penalties (Gotoh).

The original pipeline uses unit costs (mismatch 1, gap 1) and an edit distance. A
``Scoring`` object switches the pairwise step to similarity scores with separate gap
opening and extension penalties, which is what makes protein alignments biologically
meaningful.
"""

from __future__ import annotations

from dataclasses import dataclass

# BLOSUM62 as distributed by NCBI (ftp.ncbi.nlm.nih.gov/blast/matrices/BLOSUM62).
_BLOSUM62_TEXT = """
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1
"""

MATCH, DELETE, INSERT = "M", "D", "I"
_NEG_INF = float("-inf")


def _parse_matrix(text: str) -> dict[tuple[str, str], int]:
    lines = [line.split() for line in text.strip().splitlines()]
    columns = lines[0]
    matrix = {}
    for row in lines[1:]:
        for col, value in zip(columns, row[1:]):
            matrix[(row[0], col)] = int(value)
    return matrix


BLOSUM62 = _parse_matrix(_BLOSUM62_TEXT)

_PURINES = set("AG")
_PYRIMIDINES = set("CTU")


@dataclass(frozen=True)
class Scoring:
    """Similarity scoring for the pairwise step.

    ``gap_open`` is paid once per gap and ``gap_extend`` for every gap position
    (a gap of length k costs ``gap_open + k * gap_extend``); both are given as
    positive penalties.
    """

    name: str
    gap_open: float
    gap_extend: float
    matrix: dict[tuple[str, str], int] | None = None
    match: float = 0.0
    transition: float = 0.0
    transversion: float = 0.0

    @classmethod
    def blosum62(cls, gap_open: float = 10.0, gap_extend: float = 1.0) -> Scoring:
        return cls("blosum62", gap_open, gap_extend, matrix=BLOSUM62)

    @classmethod
    def dna(
        cls,
        match: float = 5.0,
        transition: float = -1.0,
        transversion: float = -4.0,
        gap_open: float = 10.0,
        gap_extend: float = 1.0,
    ) -> Scoring:
        """Nucleotides: transitions (A<->G, C<->T) are penalised less than transversions."""
        return cls("dna", gap_open, gap_extend, match=match, transition=transition, transversion=transversion)

    def score(self, a: str, b: str) -> float:
        if self.matrix is not None:
            key = (a, b)
            if key in self.matrix:
                return self.matrix[key]
            # Unknown letters score like X, gaps in a consensus like the stop column
            return self.matrix[("X" if a not in "-*" else "*", "X" if b not in "-*" else "*")]
        if a == b:
            return self.match
        pair = {a, b}
        if pair <= _PURINES or pair <= _PYRIMIDINES:
            return self.transition
        return self.transversion


def gotoh(seq1: str, seq2: str, scoring: Scoring) -> tuple[float, list[str]]:
    """Global alignment with affine gaps; returns the best score and the M/D/I path.

    Three matrices: M ends in a (mis)match, X in a gap in seq2 (DELETE, seq1 residue
    against a gap), Y in a gap in seq1 (INSERT). Ties prefer M, then DELETE, then
    INSERT, the same order as the unit-cost routine, so results stay deterministic.
    """
    n, m = len(seq1), len(seq2)
    go, ge = scoring.gap_open, scoring.gap_extend
    M = [[_NEG_INF] * (m + 1) for _ in range(n + 1)]
    X = [[_NEG_INF] * (m + 1) for _ in range(n + 1)]
    Y = [[_NEG_INF] * (m + 1) for _ in range(n + 1)]
    M[0][0] = 0.0
    for i in range(1, n + 1):
        X[i][0] = -go - i * ge
    for j in range(1, m + 1):
        Y[0][j] = -go - j * ge

    for i in range(1, n + 1):
        a = seq1[i - 1]
        for j in range(1, m + 1):
            s = scoring.score(a, seq2[j - 1])
            M[i][j] = max(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1]) + s
            X[i][j] = max(M[i - 1][j] - go - ge, X[i - 1][j] - ge, Y[i - 1][j] - go - ge)
            Y[i][j] = max(M[i][j - 1] - go - ge, Y[i][j - 1] - ge, X[i][j - 1] - go - ge)

    tables = {MATCH: M, DELETE: X, INSERT: Y}
    state = max((MATCH, DELETE, INSERT), key=lambda k: (tables[k][n][m], -"MDI".index(k)))
    best = tables[state][n][m]

    path = []
    i, j = n, m
    while i > 0 or j > 0:
        path.append(state)
        if state == MATCH:
            s = scoring.score(seq1[i - 1], seq2[j - 1])
            prev = M[i][j] - s
            i, j = i - 1, j - 1
            state = _pick(prev, ((MATCH, M[i][j]), (DELETE, X[i][j]), (INSERT, Y[i][j])))
        elif state == DELETE:
            value = X[i][j]
            i -= 1
            state = _pick(
                value, ((MATCH, M[i][j] - go - ge), (DELETE, X[i][j] - ge), (INSERT, Y[i][j] - go - ge))
            )
        else:
            value = Y[i][j]
            j -= 1
            state = _pick(
                value, ((MATCH, M[i][j] - go - ge), (INSERT, Y[i][j] - ge), (DELETE, X[i][j] - go - ge))
            )
        if i == 0 and j > 0:
            state = INSERT
        elif j == 0 and i > 0:
            state = DELETE
    path.reverse()
    return best, path


def _pick(target: float, options) -> str:
    for state, value in options:
        if abs(value - target) < 1e-9:
            return state
    return options[0][0]


def identity_distance(seq1: str, seq2: str, path: list[str]) -> float:
    """1 - identity over aligned residue pairs; 1.0 when nothing aligns."""
    i = j = same = pairs = 0
    for op in path:
        if op == MATCH:
            pairs += 1
            same += seq1[i] == seq2[j]
        if op != INSERT:
            i += 1
        if op != DELETE:
            j += 1
    return 1.0 if pairs == 0 else 1.0 - same / pairs
