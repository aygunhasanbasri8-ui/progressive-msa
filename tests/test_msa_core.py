from progressive_msa import ProgressiveMSA


def test_needleman_wunsch_distance_identical_sequences():
    msa = ProgressiveMSA({}, verbose=False)
    assert msa.needleman_wunsch_distance("ACTG", "ACTG") == 0


def test_needleman_wunsch_distance_counts_edits():
    msa = ProgressiveMSA({}, verbose=False)
    assert msa.needleman_wunsch_distance("ACTG", "ACTC") == 1


def test_calculate_distance_matrix_is_symmetric():
    sequences = {"a": "ACTG", "b": "ACTC", "c": "ACCGT"}
    msa = ProgressiveMSA(sequences, verbose=False)
    matrix = msa.calculate_distance_matrix()

    for i in range(len(matrix)):
        for j in range(len(matrix)):
            assert matrix[i][j] == matrix[j][i]
        assert matrix[i][i] == 0


def test_run_produces_aligned_sequences_of_equal_length():
    sequences = {"Seq1": "ACTG", "Seq2": "ACTC", "Seq3": "ACCGT"}
    msa = ProgressiveMSA(sequences, verbose=False)
    result = msa.run()

    assert set(result.keys()) == set(sequences.keys())
    lengths = {len(seq) for seq in result.values()}
    assert len(lengths) == 1


def test_run_preserves_ungapped_residues():
    sequences = {"Seq1": "ACTG", "Seq2": "ACTC"}
    msa = ProgressiveMSA(sequences, verbose=False)
    result = msa.run()

    for name, aligned in result.items():
        assert aligned.replace("-", "") == sequences[name]
