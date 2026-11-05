import re
from argparse import Namespace
from io import StringIO

import pytest

import batch_launcher


def read_fragments(directory, base="fragment", suffix=".in"):
    fragment_paths = sorted(directory.glob(f"{base}*{suffix}"))
    return [path.read_text() for path in fragment_paths]


def test_fasta_record_sizer_counts_sequence_bases_only():
    record = [">read-1\n", " AC \n", "GT\n"]

    assert batch_launcher.fastaRecordSizer(record) == 4


def test_fastq_record_sizer_counts_sequence_not_quality():
    record = ["@read-1\n", "ACGT\n", "+read-1\n", "!!!!\n"]

    assert batch_launcher.fastqRecordSizer(record) == 4


def test_fragmenting_fasta_preserves_records_and_content(tmp_path):
    fasta = (
        ">read-1\nAAA\n"
        ">read-2\nCCCC\n"
        ">read-3\nGGGGG\n"
    )
    file_type = batch_launcher.FragmentableFileType(
        "fasta",
        sepRE=re.compile(r"^>"),
        sizer=batch_launcher.fastaRecordSizer,
    )

    fragment_count = batch_launcher.fragmentInputStreamBySize(
        StringIO(fasta), str(tmp_path), chunk=4, fileType=file_type,
        frag_base="fragment",
    )
    fragments = read_fragments(tmp_path)

    assert fragment_count == 3
    assert fragments == [
        ">read-1\nAAA\n",
        ">read-2\nCCCC\n",
        ">read-3\nGGGGG\n",
    ]
    assert "".join(fragments) == fasta


def test_fragmenting_by_record_count_keeps_records_whole(tmp_path):
    fasta = ">read-1\nA\n>read-2\nC\n>read-3\nG\n>read-4\nT\n>read-5\nN\n"
    file_type = batch_launcher.FragmentableFileType(
        "fasta", sepRE=re.compile(r"^>")
    )

    fragment_count = batch_launcher.fragmentInputStreamBySize(
        StringIO(fasta), str(tmp_path), chunk=2, fileType=file_type,
        frag_base="fragment", splitOnSize=False,
    )
    fragments = read_fragments(tmp_path)

    assert fragment_count == 3
    assert [fragment.count(">") for fragment in fragments] == [2, 2, 1]
    assert "".join(fragments) == fasta


def test_fragmenting_single_record_writes_one_fragment(tmp_path):
    fasta = ">only-read\nACGT\n"
    file_type = batch_launcher.FASTA

    fragment_count = batch_launcher.fragmentInputStreamBySize(
        StringIO(fasta), str(tmp_path), chunk=100, fileType=file_type,
        frag_base="fragment",
    )

    assert fragment_count == 1
    assert read_fragments(tmp_path) == [fasta]


def test_get_file_type_rejects_conflicting_split_options():
    options = Namespace(
        pattern=r"^>",
        numLines=4,
        infileType=None,
    )

    with pytest.raises(Exception, match="Conflicting options"):
        batch_launcher.getFileType(options, "reads.fasta")


def test_option_hashes_map_positional_and_flagged_arguments():
    options = Namespace(
        inputFlag="2",
        outputFlags=["-o", "4"],
        threadsFlag="-t",
        prefixFlag=None,
        rel_paths=["-d"],
    )

    positional_args, flagged_args = batch_launcher.getOptionHashes(options)

    assert positional_args == {2: "in", 4: "out"}
    assert flagged_args == {
        "-o": "out",
        "-t": "threads",
        "-d": "rel",
    }


def test_local_execution_is_the_default():
    options = batch_launcher.set_up_cmd_parser().parse_args([])

    assert options.queue == batch_launcher.LOCAL
