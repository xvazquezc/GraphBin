#!/usr/bin/env python3

import csv
import logging
import os
import re
import subprocess
import sys

from cogent3.parse.fasta import MinimalFastaParser
from igraph import *

from graphbin.bidirectionalmap.bidirectionalmap import BidirectionalMap


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Anuradha Wickramarachchi", "Yu Lin"]
__license__ = "BSD-3"
__version__ = "1.7.4"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"


logger = logging.getLogger(f"GraphBin {__version__}")


def get_initial_binning_result(
    n_bins,
    bins_list,
    contig_bins_file,
    contigs_map_rev,
    graph_to_contig_map_rev,
    delimiter,
):
    logger.info("Obtaining the initial binning result")

    bins = [set() for x in range(n_bins)]
    n_skipped = 0

    try:
        with open(contig_bins_file) as contig_bins:
            readCSV = csv.reader(contig_bins, delimiter=delimiter)
            for row in readCSV:
                # Extract short contig name (strip FASTA description fields)
                contig_name = row[0].split()[0]

                if contig_name not in graph_to_contig_map_rev:
                    n_skipped += 1
                    continue

                contig_num = contigs_map_rev[
                    int(graph_to_contig_map_rev[contig_name])
                ]

                bin_num = bins_list.index(row[1])
                bins[bin_num].add(contig_num)

    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that you have provided the correct assembler type and the correct path to the binning result file in the correct format."
        )
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    if n_skipped > 0:
        logger.warning(
            f"{n_skipped} contigs in the binning result were not found in the assembly graph and were skipped."
        )

    return bins


def parse_graph(assembly_graph_file, original_contigs):
    node_count = 0

    graph_contigs = {}

    links = []

    my_map = BidirectionalMap()

    try:
        # Get links from .gfa file
        with open(assembly_graph_file) as file:
            for line in file.readlines():
                line = line.strip()

                # Identify lines with link information
                if line.startswith("L"):
                    link = []

                    strings = line.split("\t")

                    start_1 = "NODE_"
                    end_1 = "_length"

                    link1 = int(
                        re.search("%s(.*)%s" % (start_1, end_1), strings[1]).group(1)
                    )

                    start_2 = "NODE_"
                    end_2 = "_length"

                    link2 = int(
                        re.search("%s(.*)%s" % (start_2, end_2), strings[3]).group(1)
                    )

                    link.append(link1)
                    link.append(link2)
                    links.append(link)

                elif line.startswith("S"):
                    strings = line.split()

                    start = "NODE_"
                    end = "_length"

                    contig_num = int(
                        re.search("%s(.*)%s" % (start, end), strings[1]).group(1)
                    )

                    my_map[node_count] = int(contig_num)

                    graph_contigs[contig_num] = strings[2]

                    node_count += 1

        logger.info(f"Total number of contigs available: {node_count}")

        contigs_map = my_map
        contigs_map_rev = my_map.inverse

        # Create graph
        assembly_graph = Graph()

        # Add vertices
        assembly_graph.add_vertices(node_count)

        # Create list of edges
        edge_list = []

        for i in range(node_count):
            assembly_graph.vs[i]["id"] = i
            assembly_graph.vs[i]["label"] = str(contigs_map[i])

        # Iterate links
        for link in links:
            # Remove self loops
            if link[0] != link[1]:
                # Add edge to list of edges
                edge_list.append((contigs_map_rev[link[0]], contigs_map_rev[link[1]]))

        # Add edges to the graph
        assembly_graph.add_edges(edge_list)
        assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

    except BaseException as err:
        logger.error(f"Unexpected {err}")
        logger.error(
            "Please make sure that the correct path to the assembly graph file is provided."
        )
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    logger.info(f"Total number of edges in the assembly graph: {len(edge_list)}")

    # Map original contig IDs to contig IDS of assembly graph
    # --------------------------------------------------------

    graph_to_contig_map = BidirectionalMap()

    # Try matching by sequence content. Keep support for duplicate sequences by
    # storing all original IDs per sequence and consuming them one-by-one.
    seq_to_original = {}
    for name, seq in original_contigs.items():
        seq_to_original.setdefault(seq, []).append(name)

    for n, m in graph_contigs.items():
        if m in seq_to_original and seq_to_original[m]:
            graph_to_contig_map[n] = seq_to_original[m].pop()

    # Fall back to positional matching if sequence matching produced poor results
    if len(graph_to_contig_map) < len(graph_contigs) * 0.5:
        logger.warning(
            f"Sequence matching mapped only {len(graph_to_contig_map)}/{len(graph_contigs)} contigs. "
            f"Falling back to positional matching."
        )
        graph_to_contig_map = BidirectionalMap()
        original_keys = list(original_contigs.keys())
        for i, n in enumerate(graph_contigs):
            if i < len(original_keys):
                graph_to_contig_map[n] = original_keys[i]

    logger.info(f"Mapped {len(graph_to_contig_map)}/{len(graph_contigs)} contigs to original IDs")

    return assembly_graph, graph_to_contig_map, contigs_map, node_count


def write_output(
    output_path,
    prefix,
    final_bins,
    contigs_file,
    graph_to_contig_map,
    bins,
    contigs_map,
    bins_list,
    delimiter,
    node_count,
    remove_labels,
    non_isolated,
):
    logger.info("Writing the Final Binning result to file")

    output_bins = []
    graph_to_contig_map_rev = graph_to_contig_map.inverse
    contigs_map_rev = contigs_map.inverse

    output_bins_path = f"{output_path}{prefix}bins/"
    output_file = f"{output_path}{prefix}graphbin_output.csv"

    if not os.path.isdir(output_bins_path):
        subprocess.run(f"mkdir -p {output_bins_path}", shell=True)

    bin_files = {}

    for bin_name in set(final_bins.values()):
        bin_files[bin_name] = open(
            f"{output_bins_path}{prefix}bin_{bin_name}.fasta", "w+"
        )

    n_missing_fasta = 0
    n_missing_output = 0
    n_missing_unbinned = 0

    for label, seq in MinimalFastaParser(
        contigs_file, label_to_name=lambda x: x.split()[0]
    ):
        if label not in graph_to_contig_map_rev:
            n_missing_fasta += 1
            continue

        graph_contig_id = graph_to_contig_map_rev[label]

        if graph_contig_id not in contigs_map_rev:
            n_missing_fasta += 1
            continue

        contig_num = contigs_map_rev[graph_contig_id]

        if contig_num in final_bins:
            bin_files[final_bins[contig_num]].write(f">{label}\n{seq}\n")

    # Close output files
    for c in set(final_bins.values()):
        bin_files[c].close()

    for b in range(len(bins)):
        for contig in bins[b]:
            if contig not in contigs_map:
                n_missing_output += 1
                continue

            graph_contig_id = contigs_map[contig]

            if graph_contig_id not in graph_to_contig_map:
                n_missing_output += 1
                continue

            line = []
            line.append(graph_to_contig_map[graph_contig_id])
            line.append(bins_list[b])
            output_bins.append(line)

    with open(output_file, mode="w") as out_file:
        output_writer = csv.writer(
            out_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL
        )
        for row in output_bins:
            output_writer.writerow(row)

    logger.info(f"Final binning results can be found in {output_bins_path}")

    unbinned_contigs = []

    for i in range(node_count):
        if i in remove_labels or i not in non_isolated:
            if i not in contigs_map:
                n_missing_unbinned += 1
                continue

            graph_contig_id = contigs_map[i]

            if graph_contig_id not in graph_to_contig_map:
                n_missing_unbinned += 1
                continue

            line = []
            line.append(graph_to_contig_map[graph_contig_id])
            unbinned_contigs.append(line)

    if n_missing_fasta > 0:
        logger.warning(
            f"Skipped {n_missing_fasta} contigs while writing FASTA bins because they could not be mapped to original IDs."
        )

    if n_missing_output > 0:
        logger.warning(
            f"Skipped {n_missing_output} contigs while writing graphbin_output.csv because they could not be mapped to original IDs."
        )

    if n_missing_unbinned > 0:
        logger.warning(
            f"Skipped {n_missing_unbinned} contigs while writing graphbin_unbinned.csv because they could not be mapped to original IDs."
        )

    if len(unbinned_contigs) != 0:
        unbinned_file = f"{output_path}{prefix}graphbin_unbinned.csv"

        with open(unbinned_file, mode="w") as out_file:
            output_writer = csv.writer(
                out_file, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL
            )

            for row in unbinned_contigs:
                output_writer.writerow(row)

        logger.info(f"Unbinned contigs can be found at {unbinned_file}")


def get_contig_descriptors(contigs_file):
    original_contigs = {}
    contig_descriptions = {}

    for label, seq in MinimalFastaParser(contigs_file):
        name = label.split()[0]
        original_contigs[name] = seq
        contig_descriptions[name] = label

    return original_contigs
