import functools
import random

from progressive_msa import BLOSUM62, ProgressiveMSA, Scoring
from progressive_msa.cli import main
from progressive_msa.scoring import gotoh


def _brute_force(a, b, s):
    @functools.lru_cache(None)
    def best(i, j, last):
        if i == len(a) and j == len(b):
            return 0.0
        options = []
        if i < len(a) and j < len(b):
            options.append(s.score(a[i], b[j]) + best(i + 1, j + 1, "M"))
        if i < len(a):
            options.append(-(s.gap_extend + (s.gap_open if last != "D" else 0)) + best(i + 1, j, "D"))
        if j < len(b):
            options.append(-(s.gap_extend + (s.gap_open if last != "I" else 0)) + best(i, j + 1, "I"))
        return max(options)

    return best(0, 0, None)


def test_blosum62_is_the_ncbi_table():
    assert len(BLOSUM62) == 24 * 24
    assert BLOSUM62[("W", "W")] == 11 and BLOSUM62[("C", "C")] == 9 and BLOSUM62[("A", "R")] == -1
    assert all(BLOSUM62[(a, b)] == BLOSUM62[(b, a)] for a, b in BLOSUM62)


def test_gotoh_matches_brute_force():
    rnd = random.Random(3)
    for scoring, alphabet in ((Scoring.blosum62(), "ARNDCQEGHILKMFPSTWYV"), (Scoring.dna(), "ACGT")):
        for _ in range(150):
            a = "".join(rnd.choice(alphabet) for _ in range(rnd.randint(0, 8)))
            b = "".join(rnd.choice(alphabet) for _ in range(rnd.randint(0, 8)))
            score, path = gotoh(a, b, scoring)
            assert abs(score - _brute_force(a, b, scoring)) < 1e-9
            assert path.count("M") + path.count("D") == len(a)
            assert path.count("M") + path.count("I") == len(b)


def test_affine_gap_keeps_a_deletion_in_one_piece():
    full = "MKTAYIAKQRQISFVKSHFSRQ"
    seqs = {"a": full, "b": full.replace("QRQIS", "Q"), "c": "MKTAYIAKERQISFVKSHFSRQ"}

    unit = ProgressiveMSA(seqs, verbose=False).run()["b"]
    affine = ProgressiveMSA(seqs, verbose=False, scoring=Scoring.blosum62()).run()["b"]

    assert unit.count("-") == affine.count("-") == 4
    assert "----" in affine and "----" not in unit


def test_dna_transition_costs_less_than_transversion():
    s = Scoring.dna()
    assert s.score("A", "G") > s.score("A", "C")
    assert s.score("C", "T") == s.score("A", "G")


def test_cli_scoring_option(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(">a\nMKTAYIAKQRQISFVKSHFSRQ\n>b\nMKTAYIAKQFVKSHFSRQ\n")
    out = tmp_path / "out.fa"
    assert main(["align", str(src), "-o", str(out), "--scoring", "blosum62"]) == 0
    assert "MKTAYIAKQ----FVKSHFSRQ" in out.read_text()
