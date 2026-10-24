import pytest

from progressive_msa import FastaError, format_clustal, format_fasta, parse_fasta
from progressive_msa.cli import main


def test_parse_fasta_joins_lines_and_keeps_first_header_word():
    text = ">s1 Homo sapiens\nac-tg\n\n>s2\nAC\nTC\n"
    assert parse_fasta(text.splitlines()) == {"s1": "ACTG", "s2": "ACTC"}


@pytest.mark.parametrize(
    "text",
    ["ACGT\n", ">\nACGT\n", ">a\nAC\n>a\nGT\n", ">a\n>b\nGT\n", ""],
)
def test_parse_fasta_rejects_bad_input(text):
    with pytest.raises(FastaError):
        parse_fasta(text.splitlines())


def test_formatters_wrap_long_rows():
    alignment = {"a": "AC-T" * 20, "bb": "ACGT" * 20}
    fasta = format_fasta(alignment)
    assert max(len(line) for line in fasta.splitlines()) == 60

    clustal = format_clustal(alignment, block=50).splitlines()
    assert clustal[0].startswith("CLUSTAL")
    assert clustal[3] == "a     " + alignment["a"][:50]
    assert clustal[5].strip().startswith("**")


def test_cli_writes_alignment_in_input_order(tmp_path):
    src = tmp_path / "in.fa"
    src.write_text(">Seq1\nACTG\n>Seq2\nACTC\n>Seq3\nACCGT\n")
    out = tmp_path / "out.fa"

    assert main(["align", str(src), "-o", str(out)]) == 0
    assert out.read_text() == ">Seq1\nACTG-\n>Seq2\nACTC-\n>Seq3\nACCGT\n"


def test_cli_reports_missing_file(capsys):
    assert main(["align", "does-not-exist.fa"]) == 2
    assert "msa:" in capsys.readouterr().err
