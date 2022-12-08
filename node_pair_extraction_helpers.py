#!/pkg/python/3.7.4/bin/python3
from collections import defaultdict
from file_helpers import *
from ortholog_helpers import *
from graph_helpers import *

def extract_node_pairs(seeds, ignore_deg_1=False, gtag1='', gtag2=''):
    m2m_pairs = seeds_to_m2m(seeds)
    pairs = extract_node_pairs_from_m2m(m2m_pairs)

    if ignore_deg_1:
        pairs = get_rid_of_deg1_pairs(pairs, gtag1, gtag2)

    return pairs

def extract_big_patch_alignment_from_m2m(m2m_pairs, mindvs=1, minratio=1):
    node_pair_voting = create_node_pair_voting(m2m_pairs)
    node_favorite_pairs = create_node_favorite_pairs(node_pair_voting)
    
    # remove nodes with more than one favorite
    no_ties_node_favorite_pairs = dict()
    for node, favorites in node_favorite_pairs.items():
        if len(favorites) == 0:
            raise AssertionError('len of favorites should never be 0')
        elif len(favorites) == 1:
            no_ties_node_favorite_pairs[node] = favorites
    node_favorite_pairs = no_ties_node_favorite_pairs

    # remove nodes whose favorites are not selective enough
    # this is separate from removing ties. we want to remove ties no matter what, even if minratio is 0
    selective_node_favorite_pairs = dict()
    i = 0
    for node, favorites in node_favorite_pairs.items():
        assert len(favorites) == 1, 'ties weren\'t removed properly'
        voting = node_pair_voting[node].items()
        nodes_in_order = sorted(voting, key=(lambda data : data[1]), reverse=True)
        if len(nodes_in_order) == 1: # there is no ratio between fav and second fav here
            selective_node_favorite_pairs[node] = favorites
        else:
            fav_node = nodes_in_order[0]
            second_fav_node = nodes_in_order[1]
            assert fav_node[1] != second_fav_node[1], 'ties should have been removed'
            assert fav_node[1] > second_fav_node[1], 'nodes in order not correct?'
            ratio = fav_node[1] / second_fav_node[1]
            if ratio >= minratio:
                selective_node_favorite_pairs[node] = favorites
    node_favorite_pairs = selective_node_favorite_pairs
            
    # create alignment while accounting for mindvs
    alignment = []
    
    for node, favorites in node_favorite_pairs.items():
        for fav in favorites:
            if node == fav:
                exit('node equals fav')

            if node < fav: # only process in one direction to avoid duplicates
                if fav in node_favorite_pairs and node in node_favorite_pairs[fav]:
                    dvs = min(node_pair_voting[node][fav], node_pair_voting[fav][node])

                    if dvs >= mindvs:
                        alignment.append((deaug(node), deaug(fav)))

    return alignment
    
def extract_big_patch_alignment_from_seeds(seeds, mindvs=1, minratio=1):
    m2m_pairs = seeds_to_m2m(seeds)
    return extract_big_patch_alignment_from_m2m(m2m_pairs, mindvs=mindvs, minratio=minratio)

# extracts node pairs from many2many alignments (.aln files)
def extract_node_pairs_from_m2m(m2m_pairs):
    node_pair_voting = create_node_pair_voting(m2m_pairs)
    node_favorite_pairs = create_node_favorite_pairs(node_pair_voting)
    output_pairs = create_output_pairs(node_favorite_pairs)
    return output_pairs

def aug(node, n):
    return f'{n}_{node}'

def deaug(auged_node):
    return '_'.join(auged_node.split('_')[1:])

def print_output_pairs(output_pairs):
    print('\n'.join([f'{deaug(node1)} {deaug(node2)}' for node1, node2 in output_pairs]))

def get_rid_of_deg1_pairs(pairs, gtag1, gtag2):
    graph_path1 = get_graph_path(gtag1)
    nodes1 = read_in_nodes(graph_path1)
    adj_set1 = read_in_adj_set(graph_path1)
    graph_path2 = get_graph_path(gtag2)
    nodes2 = read_in_nodes(graph_path2)
    adj_set2 = read_in_adj_set(graph_path2)
    deg1_nodes1 = [node for node in nodes1 if len(adj_set1[node]) == 1]
    deg1_nodes2 = [node for node in nodes2 if len(adj_set2[node]) == 1]
    new_pairs = []

    for node1, node2 in pairs:
        if node1 not in deg1_nodes1 and node2 not in deg1_nodes2:
            new_pairs.append((node1, node2))

    return new_pairs

def create_node_pair_voting(m2m_pairs):
    def add_to_voting(node1, node2):
        if node1 not in node_pair_voting:
            node_pair_voting[node1] = defaultdict(int)

        if node2 not in node_pair_voting:
            node_pair_voting[node2] = defaultdict(int)

        node_pair_voting[node1][node2] += 1
        node_pair_voting[node2][node1] += 1

    node_pair_voting = dict()

    for s1_node, s2_node in m2m_pairs:
        add_to_voting(aug(s1_node, 1), aug(s2_node, 2))

    return node_pair_voting

def seeds_to_m2m(seeds):
    # has to be list, not set, because we want duplicates (they count towards the vote)
    m2m_pairs = list()

    for graphlet_id, s1_index, s2_index in seeds:
        for s1_node, s2_node in zip(s1_index, s2_index):
            m2m_pairs.append((s1_node, s2_node))

    return m2m_pairs

def create_node_favorite_pairs(node_pair_voting):
    node_favorite_pairs = defaultdict(set)

    for base, votes in node_pair_voting.items():
        max_count = max([count for count in votes.values()])

        for node, count in votes.items():
            if count == max_count:
                node_favorite_pairs[base].add(node)

    return node_favorite_pairs

def create_output_pairs(node_favorite_pairs):
    output_pairs = set()

    for node, favorites in node_favorite_pairs.items():
        for fav in favorites:
            if node == fav:
                exit('node equals fav')

            if node < fav: # only process in one direction to avoid duplicates
                if node in node_favorite_pairs[fav]:
                    output_pairs.add((deaug(node), deaug(fav)))

    return output_pairs

def is_well_formed_alignment(alignment):
    return all(len(pair) == 2 for pair in alignment)

if __name__ == '__main__':
    path = get_data_path('mcl/alignments/syeast0-syeast25-5000.txt')
    print(path)
    m2m_pairs = read_in_slashes_m2m(path)
    node_pairs = extract_node_pairs_from_m2m(m2m_pairs)
    orthopairs = get_orthopairs_list(node_pairs, SelfOrthos())
    print(f'{len(orthopairs)} / {len(node_pairs)}')
