import csv
import subprocess

from pathlib import Path

import pytest

from graphbin.bidirectionalmap.bidirectionalmap import BidirectionalMap
from graphbin.parsers.megahit_parser import write_output


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Gavin Huttley"]
__license__ = "BSD-3"
__version__ = "1.7.4"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Development"


TEST_ROOTDIR = Path(__file__).parent
EXEC_ROOTDIR = Path(__file__).parent.parent


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("tmp")


@pytest.fixture(autouse=True)
def workingdir(tmp_dir, monkeypatch):
    """set the working directory for all tests"""
    monkeypatch.chdir(tmp_dir)


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0

    Parameters
    ----------

    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def test_graphbin_version():
    """test graphbin version"""
    cmd = "graphbin --version"
    exec_command(cmd)


def test_graphbin_on_spades_dataset(tmp_dir):
    """test graphbin on spades assembly"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_metaSPAdes"
    graph = dir_name / "assembly_graph_with_scaffolds.gfa"
    contigs = dir_name / "contigs.fasta"
    paths = dir_name / "contigs.paths"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler spades --graph {graph} --contigs {contigs} --paths {paths} --binned {binned} --output {tmp_dir}"
    exec_command(cmd)


def test_graphbin_on_sga_dataset(tmp_dir):
    """test graphbin on sga assembly"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_SGA"
    graph = dir_name / "default-graph.asqg"
    contigs = dir_name / "default-contigs.fa"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler sga --graph {graph} --contigs {contigs} --binned {binned} --output {tmp_dir}"
    exec_command(cmd)


def test_graphbin_on_megahit_dataset(tmp_dir):
    """test graphbin on megahit assembly"""
    dir_name = TEST_ROOTDIR / "data" / "ESC_MEGAHIT"
    graph = dir_name / "final.gfa"
    contigs = dir_name / "final.contigs.fa"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler megahit --graph {graph} --contigs {contigs} --binned {binned} --output {tmp_dir}"
    exec_command(cmd)


def test_megahit_write_output_skips_unmapped_contigs(tmp_dir):
    """MEGAHIT output writer should not crash when some contigs are unmapped."""
    contigs_file = tmp_dir / "contigs.fasta"
    contigs_file.write_text(
        ">c0\nAAAA\n>c1\nCCCC\n>c2\nGGGG\n",
        encoding="utf-8",
    )

    graph_to_contig_map = BidirectionalMap()
    graph_to_contig_map[100] = "c0"
    graph_to_contig_map[101] = "c1"
    # Intentionally do not map graph contig id 102.

    contigs_map = BidirectionalMap()
    contigs_map[0] = 100
    contigs_map[1] = 101
    contigs_map[2] = 102

    write_output(
        output_path=f"{tmp_dir}/",
        prefix="",
        final_bins={0: "1", 1: "1", 2: "2"},
        contigs_file=contigs_file,
        graph_to_contig_map=graph_to_contig_map,
        bins=[[0, 1, 2]],
        contigs_map=contigs_map,
        bins_list=["1"],
        delimiter=",",
        node_count=3,
        remove_labels={2},
        non_isolated={0, 1},
    )

    output_file = tmp_dir / "graphbin_output.csv"
    assert output_file.exists()

    with output_file.open(encoding="utf-8") as handle:
        rows = list(csv.reader(handle))

    assert rows == [["c0", "1"], ["c1", "1"]]


def test_graphbin_on_flye_dataset(tmp_dir):
    """test graphbin on flye assembly"""
    dir_name = TEST_ROOTDIR / "data" / "1Y3B_Flye"
    graph = dir_name / "assembly_graph.gfa"
    contigs = dir_name / "assembly.fasta"
    paths = dir_name / "assembly_info.txt"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler flye --graph {graph} --contigs {contigs} --paths {paths} --binned {binned} --output {tmp_dir}"
    exec_command(cmd)


def test_graphbin_on_canu_dataset(tmp_dir):
    """test graphbin on canu assembly"""
    dir_name = TEST_ROOTDIR / "data" / "1Y3B_Canu"
    graph = dir_name / "1y3b.contigs.gfa"
    contigs = dir_name / "1y3b.contigs.fasta"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler canu --graph {graph} --contigs {contigs} --binned {binned} --output {tmp_dir}"
    exec_command(cmd)


def test_graphbin_on_miniasm_dataset(tmp_dir):
    """test graphbin on miniasm assembly"""
    dir_name = TEST_ROOTDIR / "data" / "1Y3B_Miniasm"
    graph = dir_name / "reads.gfa"
    contigs = dir_name / "unitigs.fasta"
    binned = dir_name / "initial_binning_res.csv"
    cmd = f"graphbin --assembler miniasm --graph {graph} --contigs {contigs} --binned {binned} --output {tmp_dir}"
    exec_command(cmd)
