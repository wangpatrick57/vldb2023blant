#!/pkg/python/3.7.4/bin/python3
import re
import math
import sys
from file_helpers import *
from graph_helpers import *

NODE_LIMIT=20_000
EDGE_LIMIT=400_000

def clean_tel(graph_path):
    from graph_helpers import get_canon_edge
    
    tel = []
    
    with open(graph_path, 'r') as graph_file:
        added_canon_edges = set()

        for i, line in enumerate(graph_file):
            node1, node2, time = re.split('[\s\t]', line.strip())
            time = int(time)

            if node1 == node2:
                continue

            if get_canon_edge(node1, node2) in added_canon_edges:
                continue
            
            tel.append((node1, node2, time))
            added_canon_edges.add(get_canon_edge(node1, node2))

    with open(graph_path, 'w') as graph_file:
        write_to_file(tel_to_str(tel), graph_path)

def read_in_temporal_el(graph_path, edge_limit=None):
    from graph_helpers import get_canon_edge
    
    with open(graph_path, 'r') as graph_file:
        tel = []
        num_edges = 0
        added_canon_edges = set()

        for i, line in enumerate(graph_file):
            node1, node2, time = re.split('[\s\t]', line.strip())
            time = int(time)
            
            assert node1 != node2
            assert get_canon_edge(node1, node2) not in added_canon_edges

            added_canon_edges.add(get_canon_edge(node1, node2))
            tel.append((node1, node2, time))
            num_edges += 1

            if edge_limit != None and num_edges >= edge_limit:
                break

        return tel

def get_tel_sections(tel, length_ratios, offset_ratios):
    assert len(length_ratios) == len(offset_ratios), 'length_ratios must be the same len as offset_ratios' 
    times = [time for node1, node2, time in tel]
    min_time = min(times)
    max_time = max(times)
    total_length = max_time - min_time
    els = []

    for length_ratio, offset_ratio in zip(length_ratios, offset_ratios):
        length = length_ratio * total_length
        offset = offset_ratio * total_length
        start_time = min_time + offset
        end_time = start_time + length
        el = get_el_in_interval(tel, start_time, end_time)
        els.append(el)

    return els

def get_el_in_interval(tel, start_time, end_time):
    el = []

    for node1, node2, time in tel:
        if start_time <= time < end_time:
            el.append((node1, node2))

    return el

# start time is inclusive, end time is exclusive
def read_in_xel_in_interval(graph_path, start_time, end_time, is_tel=True):
    with open(graph_path, 'r') as graph_file:
        xel = []

        for line in graph_file:
            node1, node2, time = re.split('[\s\t]', line.strip())
            time = int(time)

            if time >= end_time:
                break

            if time >= start_time:
                edge = [source, target]

                if is_tel:
                    edge.append(str(int(float(time))))

                edge = tuple(edge)
                xel.append(edge)

        return xel

def nodes_of_tel(tel):
    nodes = set()
    
    for node1, node2, time in tel:
        nodes.add(node1)
        nodes.add(node2)

    return nodes

def tel_to_str(tel):
    return '\n'.join([f'{node1}\t{node2}\t{time}' for node1, node2, time in tel])

def get_el_from_tel_with_limits(tel, num_skipped_edges, node_limit=NODE_LIMIT, edge_limit=EDGE_LIMIT, density_limit=20, total_edge_ratio_limit=10/11):
    from graph_helpers import nodes_of_el

    tel_max_edges = len(tel) * total_edge_ratio_limit
    max_edges = min(edge_limit, tel_max_edges)
    num_edges = 0
    nodes = set()
    total_nodes = len(nodes_of_tel(tel))
    el = []
    tel_slice = tel[num_skipped_edges:]

    for node1, node2, time in tel_slice:
        if num_edges >= max_edges:
            break

        num_nodes = len(nodes)
        
        if num_nodes >= node_limit:
            break

        if num_nodes >= total_nodes / 100:
            if num_edges / num_nodes >= density_limit:
                break

        # we don't add self loops just cuz blant doesn't allow them. however, we still count them as edges in num_edges
        if node1 != node2:
            el.append((node1, node2))
            nodes.add(node1)
            nodes.add(node2)
            
        num_edges += 1 # we intentionally count duplicates as different edges to avoid the situation where a small graph with lots of duplicates never hits the node, edge, or density limits and generates graph slices that just span the entire graph

    return el

def get_std_percents():
    return [0, 1, 3, 5]

def get_gtag_from_tgtag(tgtag, percent):
    return f'{tgtag}_s{percent}'

def gen_std_tgtag_els(tgtag):
    from graph_helpers import el_to_str
    from file_helpers import write_to_file

    tpath = get_tgraph_path(tgtag)
    tel = read_in_temporal_el(tpath, EDGE_LIMIT * 1.1)
    num_edges = len(tel)

    for percent in get_std_percents():
        num_skipped_edges = num_edges * percent // 100
        el = get_el_from_tel_with_limits(tel, num_skipped_edges)
        new_gtag = get_gtag_from_tgtag(tgtag, percent)
        new_path = get_graph_path(new_gtag)
        el_str = el_to_str(el)
        write_to_file(el_str, new_path)

def map_density_over_time(tel, granularity=100):
    times = [int(row[2]) for row in tel]
    start_time = min(times)
    end_time = max(times) + 1
    tbins = [0] * granularity

    for time in times:
        tbin_num = math.floor((time - start_time) * granularity / (end_time - start_time))
        tbins[tbin_num] += 1

    for i, lines_in_this_range in enumerate(tbins):
        this_start_time = start_time + (end_time - start_time) * i / granularity
        print(this_start_time, lines_in_this_range)

def get_tgraph_path(tgtag):
    from graph_helpers import get_base_graph_path
    return get_base_graph_path(f'snap/{tgtag}', ext='tel')

if __name__ == '__main__':
    from graph_helpers import get_paper_tprl_snap

    for tgtag in get_paper_tprl_snap():
        print(f'starting {tgtag}')
        gen_std_tgtag_els(tgtag)
        print(f'done with {tgtag}')
