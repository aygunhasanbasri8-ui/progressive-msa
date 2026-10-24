import pytest

from progressive_msa import ProgressiveMSA, Scoring, sum_of_pairs, total_column
from progressive_msa.cli import main

REFERENCE = {"a": "AC-GT", "b": "ACCGT", "c": "A--GT"}


def test_identical_alignment_scores_one():
    assert sum_of_pairs(REFERENCE, REFERENCE) == 1.0
    assert total_column(REFERENCE, REFERENCE) == 1.0


def test_shifted_column_loses_pairs_and_columns():
    test = {"a": "ACG-T", "b": "ACCGT", "c": "A--GT"}
    assert 0 < sum_of_pairs(test, REFERENCE) < 1
    assert total_column(test, REFERENCE) < 1


def test_different_sequences_are_rejected():
    with pytest.raises(ValueError):
        sum_of_pairs({"a": "AC", "b": "AC"}, {"a": "AC", "b": "AG"})


def test_neighbor_joining_recovers_additive_tree():
    # Additive distances of ((a,b),(c,d)); UPGMA would get this right too, NJ must as well
    names = ["a", "b", "c", "d"]
    msa = ProgressiveMSA(dict.fromkeys(names, "A"), verbose=False, tree="nj")
    msa.distance_matrix = [[0, 3, 8, 9], [3, 0, 9, 10], [8, 9, 0, 5], [9, 10, 5, 0]]
    tree = msa.build_guide_tree()
    groups = {frozenset(_leaves(t)) for t in tree}
    assert groups == {frozenset({"a", "b"}), frozenset({"c", "d"})}

    with pytest.raises(ValueError):
        ProgressiveMSA({}, tree="unknown")


def test_profile_mode_keeps_every_residue(tmp_path):
    seqs = {"a": "MKTAYIAKQRQ", "b": "MKTAYIAKQ", "c": "MRTAYLAKQRQ", "d": "MKSAYIAEQRQ"}
    aligned = ProgressiveMSA(seqs, verbose=False, scoring=Scoring.blosum62(), tree="nj").run()
    assert {k: v.replace("-", "") for k, v in aligned.items()} == seqs

    test = tmp_path / "t.fa"
    test.write_text("".join(f">{k}\n{v}\n" for k, v in aligned.items()))
    assert main(["score", str(test), str(test)]) == 0


def _leaves(node):
    return [node] if isinstance(node, str) else [leaf for child in node for leaf in _leaves(child)]
