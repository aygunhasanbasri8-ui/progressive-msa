"""Progressive multiple sequence alignment (Clustal Omega-style).

Implements the classic three-stage pipeline:
1. Pairwise distance matrix via Needleman-Wunsch global alignment.
2. Guide tree construction via UPGMA clustering.
3. Progressive alignment of sequences/profiles following the guide tree.
"""

from typing import Dict, List, Tuple, Union

AlignmentTree = Union[str, Tuple["AlignmentTree", "AlignmentTree"]]

MISMATCH_COST = 1
GAP_PENALTY = 1

# Traceback operations: M consumes a residue from both sides,
# D only from the first sequence, I only from the second.
MATCH, DELETE, INSERT = "M", "D", "I"


def _needleman_wunsch(seq1: str, seq2: str) -> Tuple[int, List[str]]:
    """Returns the edit distance and the alignment path between two sequences."""
    n, m = len(seq1), len(seq2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]
    tb = [[MATCH] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = i * GAP_PENALTY
        tb[i][0] = DELETE
    for j in range(1, m + 1):
        dp[0][j] = j * GAP_PENALTY
        tb[0][j] = INSERT

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            cost = 0 if seq1[i - 1] == seq2[j - 1] else MISMATCH_COST
            diag = dp[i - 1][j - 1] + cost
            up = dp[i - 1][j] + GAP_PENALTY
            left = dp[i][j - 1] + GAP_PENALTY

            # Ties go diag > up > left; changing the order changes the output alignment.
            best = min(diag, up, left)
            dp[i][j] = best
            if best == diag:
                tb[i][j] = MATCH
            elif best == up:
                tb[i][j] = DELETE
            else:
                tb[i][j] = INSERT

    path = []
    i, j = n, m
    while i > 0 or j > 0:
        op = tb[i][j]
        path.append(op)
        if op != INSERT:
            i -= 1
        if op != DELETE:
            j -= 1

    path.reverse()
    return dp[n][m], path


def _apply_path(seq: str, path: List[str], gap_op: str) -> str:
    """Inserts a gap wherever the path holds gap_op, consumes a residue otherwise."""
    residues = iter(seq)
    return "".join("-" if op == gap_op else next(residues) for op in path)


class ProgressiveMSA:
    """Aligns a set of sequences using progressive multiple sequence alignment."""

    def __init__(self, sequences: Dict[str, str], verbose: bool = True):
        """
        Args:
            sequences: Mapping of sequence name to sequence string.
            verbose: Whether to print progress information during alignment.
        """
        self.sequences = sequences
        self.names = list(sequences.keys())
        self.num_seqs = len(sequences)
        self.distance_matrix: List[List[int]] = []
        self.guide_tree: AlignmentTree = None
        self.verbose = verbose

    def _log(self, message: str) -> None:
        if self.verbose:
            print(message)

    def calculate_distance_matrix(self) -> List[List[int]]:
        """Computes the pairwise edit-distance matrix using Needleman-Wunsch."""
        self._log("1. Calculating distance matrix (Needleman-Wunsch)...")
        matrix = [[0] * self.num_seqs for _ in range(self.num_seqs)]

        for i in range(self.num_seqs):
            for j in range(i + 1, self.num_seqs):
                distance = self.needleman_wunsch_distance(
                    self.sequences[self.names[i]], self.sequences[self.names[j]]
                )
                matrix[i][j] = matrix[j][i] = distance

        self.distance_matrix = matrix
        return self.distance_matrix

    def needleman_wunsch_distance(self, seq1: str, seq2: str) -> int:
        """Returns the global-alignment edit distance between two sequences."""
        distance, _ = _needleman_wunsch(seq1, seq2)
        return distance

    def build_guide_tree(self) -> AlignmentTree:
        """Builds a guide tree over the sequences using UPGMA clustering."""
        self._log("2. Building guide tree (UPGMA)...")
        clusters: Dict[int, AlignmentTree] = dict(enumerate(self.names))
        distances: Dict[Tuple[int, int], float] = {
            (i, j): self.distance_matrix[i][j]
            for i in range(self.num_seqs)
            for j in range(self.num_seqs)
            if i != j
        }

        active = list(clusters)
        next_id = self.num_seqs

        while len(active) > 1:
            # Strict < keeps the first closest pair in scan order, which keeps the tree deterministic.
            min_dist = float("inf")
            c1 = c2 = -1
            for i, id1 in enumerate(active):
                for id2 in active[i + 1:]:
                    if distances[(id1, id2)] < min_dist:
                        min_dist, c1, c2 = distances[(id1, id2)], id1, id2

            clusters[next_id] = (clusters[c1], clusters[c2])
            for other in active:
                if other not in (c1, c2):
                    merged = (distances[(c1, other)] + distances[(c2, other)]) / 2.0
                    distances[(next_id, other)] = distances[(other, next_id)] = merged

            active.remove(c1)
            active.remove(c2)
            active.append(next_id)
            next_id += 1

        self.guide_tree = clusters[active[0]]
        self._log(f"   Tree hierarchy: {self.guide_tree}")
        return self.guide_tree

    def progressive_alignment(self) -> Dict[str, str]:
        """Aligns all sequences by traversing the guide tree bottom-up."""
        self._log("3. Running progressive alignment...\n")
        return self._traverse_and_align(self.guide_tree)

    def _traverse_and_align(self, node: AlignmentTree) -> Dict[str, str]:
        if isinstance(node, str):
            return {node: self.sequences[node]}

        left, right = node
        return self._align_profiles(self._traverse_and_align(left), self._traverse_and_align(right))

    def _get_consensus(self, profile: Dict[str, str]) -> str:
        consensus = []
        for column in zip(*profile.values()):
            residues = [c for c in column if c != "-"]
            consensus.append(max(set(residues), key=residues.count) if residues else "-")
        return "".join(consensus)

    def _align_profiles(self, prof1: Dict[str, str], prof2: Dict[str, str]) -> Dict[str, str]:
        # Profiles are aligned through their consensus strings, not column by column.
        _, path = _needleman_wunsch(self._get_consensus(prof1), self._get_consensus(prof2))

        aligned = {name: _apply_path(seq, path, INSERT) for name, seq in prof1.items()}
        aligned.update({name: _apply_path(seq, path, DELETE) for name, seq in prof2.items()})
        return aligned

    def run(self) -> Dict[str, str]:
        """Runs the full pipeline and returns the aligned sequences."""
        self.calculate_distance_matrix()

        self._log("\nDistance matrix:")
        for row in self.distance_matrix:
            self._log(str(row))
        self._log("-" * 30)

        self.build_guide_tree()
        aligned_sequences = self.progressive_alignment()

        self._log("Alignment complete!")
        return aligned_sequences
