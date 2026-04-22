#!/usr/bin/env python3

import logging
import sys

from collections import deque
from concurrent.futures import ThreadPoolExecutor

from graphbin.labelpropagation.labelprop import LabelProp


__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019-2022, GraphBin Project"
__credits__ = ["Vijini Mallawaarachchi", "Anuradha Wickramarachchi", "Yu Lin"]
__license__ = "BSD-3"
__version__ = "1.7.4"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"
__status__ = "Production"

logger = logging.getLogger(f"GraphBin {__version__}")

MIN_BIN_COUNT = 10


def getClosestLabelledVertices(graph, node, binned_contigs):
    # Remove labels of ambiguous vertices
    # -------------------------------------

    queu_l = deque([graph.neighbors(node, mode="ALL")])
    visited_l = {node}
    labelled = []

    while queu_l:
        active_level = queu_l.popleft()
        is_finish = False
        visited_l.update(active_level)

        for n in active_level:
            if n in binned_contigs:
                is_finish = True
                labelled.append(n)
        if is_finish:
            return labelled
        else:
            temp = set()
            for n in active_level:
                temp.update(graph.neighbors(n, mode="ALL"))
            temp2 = [n for n in temp if n not in visited_l]
            if temp2:
                queu_l.append(temp2)
    return labelled


def graphbin_main(
    n_bins, bins, bins_list, assembly_graph, node_count, diff_threshold, max_iteration,
    nthreads=1,
):
    logger.info("Determining ambiguous vertices")

    # Build reverse map from contig to its bin for O(1) lookups
    contig_bin_map = {contig: b for b in range(n_bins) for contig in bins[b]}

    remove_by_bin = {}

    remove_labels = set()

    neighbours_have_same_label_list = set()

    for b in range(n_bins):
        for i in bins[b]:
            my_bin = b

            # Get set of closest labelled vertices with distance = 1
            closest_neighbours = assembly_graph.neighbors(i, mode="all")

            # Determine whether all the closest labelled vertices have the same label as its own
            neighbours_have_same_label = True

            neighbours_binned = False

            for neighbour in closest_neighbours:
                k = contig_bin_map.get(neighbour, -1)
                if k != -1:
                    neighbours_binned = True
                    if k != my_bin:
                        neighbours_have_same_label = False
                        break

            if not neighbours_have_same_label:
                if my_bin in remove_by_bin:
                    if len(bins[my_bin]) - len(remove_by_bin[my_bin]) >= MIN_BIN_COUNT:
                        remove_labels.add(i)
                        remove_by_bin[my_bin].append(i)
                else:
                    if len(bins[my_bin]) >= MIN_BIN_COUNT:
                        remove_labels.add(i)
                        remove_by_bin[my_bin] = [i]

            elif neighbours_binned:
                neighbours_have_same_label_list.add(i)

    for i in remove_labels:
        n = contig_bin_map.pop(i, None)
        if n is not None:
            bins[n].discard(i)

    # Further remove labels of ambiguous vertices
    binned_contigs = set().union(*bins)

    # Collect vertices needing BFS (not already confirmed same-label in pass 1)
    _pass2_tasks = [
        (b, i)
        for b in range(n_bins)
        for i in bins[b]
        if i not in neighbours_have_same_label_list
    ]

    def _bfs(b_i):
        b, i = b_i
        return b, i, getClosestLabelledVertices(assembly_graph, i, binned_contigs)

    if nthreads > 1 and _pass2_tasks:
        with ThreadPoolExecutor(max_workers=nthreads) as _executor:
            _pass2_results = list(_executor.map(_bfs, _pass2_tasks))
    else:
        _pass2_results = [_bfs(t) for t in _pass2_tasks]

    for b, i, closest_neighbours in _pass2_results:
        if closest_neighbours:
            my_bin = b
            neighbours_have_same_label = True

            for neighbour in closest_neighbours:
                k = contig_bin_map.get(neighbour, -1)
                if k != -1 and k != my_bin:
                    neighbours_have_same_label = False
                    break

            if not neighbours_have_same_label and i not in remove_labels:
                if my_bin in remove_by_bin:
                    if (
                        len(bins[my_bin]) - len(remove_by_bin[my_bin])
                        >= MIN_BIN_COUNT
                    ):
                        remove_labels.add(i)
                        remove_by_bin[my_bin].append(i)
                else:
                    if len(bins[my_bin]) >= MIN_BIN_COUNT:
                        remove_labels.add(i)
                        remove_by_bin[my_bin] = [i]

    logger.info("Removing labels of ambiguous vertices")

    # Remove labels of ambiguous vertices
    for i in remove_labels:
        n = contig_bin_map.pop(i, None)
        if n is not None:
            bins[n].discard(i)

    logger.info("Obtaining the refined binning result")

    # Get vertices which are not isolated and not in components without any labels
    # -----------------------------------------------------------------------------

    logger.info(
        "Deteremining vertices which are not isolated and not in components without any labels"
    )

    non_isolated = set()

    for component in assembly_graph.clusters():
        if any(v in binned_contigs for v in component):
            non_isolated.update(component)

    logger.info("Number of non-isolated contigs: " + str(len(non_isolated)))

    # Run label propagation
    # -----------------------

    data = []

    for contig in range(node_count):
        # Consider vertices that are not isolated

        if contig in non_isolated:
            line = []
            line.append(contig)

            b = contig_bin_map.get(contig, -1)
            if b != -1:
                line.append(b + 1)
            else:
                line.append(0)

            neighbours = assembly_graph.neighbors(contig, mode="all")

            neighs = []

            for neighbour in neighbours:
                n = []
                n.append(neighbour)
                n.append(1.0)
                neighs.append(n)

            line.append(neighs)

            data.append(line)

    # Check if initial binning result consists of contigs belonging to multiple bins

    multiple_bins = False

    for item in data:
        if type(item[1]) is int and type(item[2]) is int:
            multiple_bins = True
            break

    if multiple_bins:
        logger.error(
            "Initial binning result consists of contigs belonging to multiple bins. Please make sure that each contig in the initial binning result belongs to only one bin."
        )
        logger.info("Exiting GraphBin... Bye...!")
        sys.exit(1)

    # Label propagation

    lp = LabelProp()

    lp.load_data_from_mem(data)

    logger.info(
        "Starting label propagation with eps="
        + str(diff_threshold)
        + " and max_iteration="
        + str(max_iteration)
    )

    ans = lp.run(diff_threshold, max_iteration, show_log=True, clean_result=False)

    logger.info("Obtaining Label Propagation result")

    for l in ans:
        b = l[1] - 1
        if 0 <= b < n_bins and l[0] not in bins[b]:
            bins[b].add(l[0])
            contig_bin_map[l[0]] = b

    # Remove labels of ambiguous vertices
    # -------------------------------------

    logger.info("Determining ambiguous vertices")

    remove_by_bin = {}

    remove_labels = set()

    for b in range(n_bins):
        for i in bins[b]:
            my_bin = b

            closest_neighbours = assembly_graph.neighbors(i, mode="all")

            # Determine whether all the closest labelled vertices have the same label as its own
            neighbours_have_same_label = True

            for neighbour in closest_neighbours:
                k = contig_bin_map.get(neighbour, -1)
                if k != -1 and k != my_bin:
                    neighbours_have_same_label = False
                    break

            if not neighbours_have_same_label:
                if my_bin in remove_by_bin:
                    if len(bins[my_bin]) - len(remove_by_bin[my_bin]) >= MIN_BIN_COUNT:
                        remove_labels.add(i)
                        remove_by_bin[my_bin].append(i)
                else:
                    if len(bins[my_bin]) >= MIN_BIN_COUNT:
                        remove_labels.add(i)
                        remove_by_bin[my_bin] = [i]

    logger.info("Removing labels of ambiguous vertices")

    # Remove labels of ambiguous vertices
    for i in remove_labels:
        n = contig_bin_map.pop(i, None)
        if n is not None:
            bins[n].discard(i)

    logger.info("Obtaining the Final Refined Binning result")

    final_bins = {}

    for i in range(n_bins):
        for contig in bins[i]:
            final_bins[contig] = bins_list[i]

    return final_bins, remove_labels, non_isolated
