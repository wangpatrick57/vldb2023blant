#!/bin/python3
from all_helpers import *

def full_run(g1_index_path, g1_graph_path, g2_index_path, g2_graph_path):
    k = 8
    g1_index = get_patched_index(k, g1_index_path, g1_graph_path)
    g2_index = get_patched_index(k, g2_index_path, g2_graph_path)
    seeds = find_seeds(g1_index, g2_index, print_progress=False)
    print(seeds_to_str(seeds))

if __name__ == '__main__':
    g1_index_path = sys.argv[1]
    g1_graph_path = sys.argv[2]
    g2_index_path = sys.argv[3]
    g2_graph_path = sys.argv[4]
    full_run(g1_index_path, g1_graph_path, g2_index_path, g2_graph_path)
