#!/bin/python3
from all_helpers import *

def full_run(alignments_path, g1_graph_path, g2_graph_path):
    auto_k = (10000, 0.01)
    adj_set1 = read_in_adj_set(g1_graph_path)
    adj_set2 = read_in_adj_set(g2_graph_path)
    seeds = read_in_seeds(alignments_path)
    blocks = seeds_to_blocks(seeds)
    sagrow = SimAnnealGrow(blocks, adj_set1, adj_set2, p_func=ALWAYS_P_FUNC, s3_threshold=1)
    alignment = sagrow.run(auto_k=auto_k, silent=False)
    alignment = get_clean_alignment(alignment, adj_set1, adj_set2)
    print(alignment_to_str(alignment))

if __name__ == '__main__':
    alignments_path = sys.argv[1]
    g1_graph_path = sys.argv[2]
    g2_graph_path = sys.argv[3]
    full_run(alignments_path, g1_graph_path, g2_graph_path)

