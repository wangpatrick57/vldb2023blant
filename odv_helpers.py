#!/pkg/python/3.7.4/bin/python3
import sys
import math
import heapq
from collections import defaultdict, namedtuple
from statistics import mean
from bash_helpers import *
from graph_helpers import *
from file_helpers import *
from general_helpers import *

# WEIRD ORCA BEHAVIOR
# it seems to output all the nodes with nothing and then all the nodes with the right ODV values?

CACHED_CUM_ORBIT_COUNTS = dict()
CACHED_CUM_ORBIT_COUNTS[6] = [1, 2, 2, 2, 3, 4, 3, 3, 4, 3, 4, 4, 4, 4, 3, 4, 6, 5, 4, 5, 6, 6, 4, 4, 5, 5, 8, 4, 6, 6, 7, 5, 6, 6, 6, 5, 6, 7, 7, 5, 7, 7, 7, 6, 5, 5, 6, 8, 8, 6, 6, 8, 6, 9, 6, 6, 4, 6, 6, 8, 9, 6, 6, 8, 8, 6, 7, 7, 8, 5, 6, 6, 4, 5, 5, 7, 5, 8, 8, 7, 8, 8, 7, 9, 7, 5, 8, 8, 9, 9, 7, 8, 12, 12, 8, 10, 8, 10, 8, 10, 10, 10, 7, 9, 11, 8, 9, 13, 7, 10, 9, 10, 7, 10, 10, 10, 8, 8, 8, 8, 7, 8, 10, 9, 9, 8, 12, 12, 7, 9, 9, 6, 6, 10, 8, 8, 6, 10, 10, 11, 6, 10, 8, 6, 10, 5, 8, 8, 8, 11, 12, 7, 6, 8, 11, 10, 12, 9, 8, 11, 11, 14, 14, 7, 11, 10, 10, 10, 11, 13, 12, 14, 15, 7, 13, 14, 10, 7, 10, 12, 7, 8, 12, 5, 8, 8, 10, 8, 7, 8, 10, 11, 12, 9, 13, 10, 14, 14, 14, 10, 13, 12, 14, 13, 16, 9, 11, 14, 12, 14, 8, 11, 12, 12, 12, 11, 6, 10, 8, 10, 11, 6, 10, 10, 12, 11, 10, 8, 12, 11, 10, 10, 10, 10, 7, 11, 11, 10, 7, 12, 12, 11, 11, 12, 11, 13, 15, 14, 11, 14, 13, 14, 14, 13, 14, 14, 11, 10, 12, 12, 10, 11, 11, 7, 12, 13, 12, 9, 10, 11, 15, 15, 6, 13, 9, 6, 8, 11, 12, 10, 9, 12, 13, 14, 13, 16, 8, 10, 10, 14, 14, 7, 13, 12, 11, 9, 15, 12, 11, 15, 7, 11, 11, 11, 9, 9, 12, 9, 14, 12, 11, 16, 12, 13, 11, 13, 14, 16, 13, 14, 16, 14, 13, 16, 10, 11, 13, 14, 13, 15, 12, 14, 14, 16, 17, 15, 17, 14, 11, 16, 11, 10, 12, 10, 14, 10, 14, 13, 16, 13, 13, 15, 13, 16, 15, 15, 11, 10, 14, 11, 14, 10, 12, 12, 12, 9, 12, 9, 12, 12, 12, 8, 9, 9, 5, 9, 10, 8, 7, 11, 13, 12, 7, 9, 12, 12, 5, 8, 8, 8, 10, 10, 8, 12, 11, 8, 9, 8, 9, 16, 13, 16, 15, 11, 12, 13, 11, 16, 15, 14, 14, 12, 10, 14, 12, 14, 15, 14, 11, 10, 15, 11, 15, 13, 17, 11, 17, 12, 16, 12, 12, 16, 15, 15, 12, 13, 15, 15, 11, 8, 12, 12, 8, 10, 10, 12, 10, 11, 9, 14, 7, 14, 17, 17, 16, 12, 17, 17, 10, 15, 14, 12, 13, 14, 10, 12, 15, 12, 14, 14, 10, 9, 10, 10, 10, 12, 8, 8, 8, 5]
CACHED_CUM_ORBIT_COUNTS[5] = [1, 2, 2, 2, 3, 4, 3, 3, 4, 3, 4, 4, 4, 4, 3, 4, 6, 5, 4, 5, 6, 6, 4, 4, 5, 5, 8, 4, 6, 6, 7, 5, 6, 6, 6, 5, 6, 7, 7, 5, 7, 7, 7, 6, 5, 5, 6, 8, 8, 6, 6, 8, 6, 9, 6, 6, 4, 6, 6, 8, 9, 6, 6, 8, 8, 6, 7, 7, 8, 5, 6, 6, 4]
CACHED_CUM_ORBIT_COUNTS[4] = [1, 2, 2, 2, 3, 4, 3, 3, 4, 3, 4, 4, 4, 4, 3]

def get_odv_path(gtag, k):
    return get_data_path(f'odv/{gtag}-k{k}.odv')

def get_blantspl_path_nstr(gtag, k, nstr):
    return get_data_path(f'odv/{gtag}-k{k}-n{nstr}.splodv')

def get_blantspl_path(gtag, k, n):
    return get_data_path(f'odv/{gtag}-k{k}-n{get_abbr_num_str(n)}.splodv')

def get_cbodv_path(gtag, k, nstr):
    return get_data_path(f'odv/{gtag}-k{k}-n{nstr}.cbodv')

def gtag_to_k(gtag, override_k=None):
    from graph_helpers import is_syeast

    if override_k != None:
        return override_k
    
    if is_syeast(gtag):
        return 5
    else:
        return 4

def get_gtag_to_n_cache_path():
    return get_data_path('caches/gtag_to_n_cache.txt')
    
def read_gtag_to_n_cache():
    g2n_cache = dict()

    with open(get_gtag_to_n_cache_path(), 'r') as f:
        for line in f:
            gtag, n = line.strip().split()
            n = int(n)
            g2n_cache[gtag] = n

    return g2n_cache

def write_gtag_to_n_cache(cache):
    cache_str = '\n'.join([f'{gtag} {n}' for gtag, n in cache.items()])
    write_to_file(cache_str, get_gtag_to_n_cache_path())

def gtag_to_n(gtag):
    from graph_helpers import read_in_nodes, get_graph_path
    
    g2n_cache = read_gtag_to_n_cache()

    if gtag in g2n_cache:
        return g2n_cache[gtag]
    else:    
        nodes = read_in_nodes(get_graph_path(gtag))
        g2n_cache[gtag] = len(nodes)
        write_gtag_to_n_cache(g2n_cache)
        return len(nodes)

def two_gtags_to_k(gtag1, gtag2, override_k=None):
    assert gtag_to_k(gtag1, override_k=override_k) == gtag_to_k(gtag2, override_k=override_k)
    k = gtag_to_k(gtag1, override_k=override_k)
    return k

def two_gtags_to_n(gtag1, gtag2):
    return min(gtag_to_n(gtag1), gtag_to_n(gtag2))

def get_num_graphlets(k):
    if k == 8:
        return 11117
    elif k == 7:
        return 853
    elif k == 6:
        return 112
    elif k == 5:
        return 21
    elif k == 4:
        return 6
    elif k == 3:
        return 2
    elif k == 2:
        return 1
    else:
        return None

def get_num_graphlets_cum(k):
    if k < 2:
        return None
    elif k == 2:
        return get_num_graphlets(k)
    else:
        return get_num_graphlets(k) + get_num_graphlets_cum(k - 1)

def get_num_orbits(k):
    if k == 8:
        return 72489
    elif k == 7:
        return 4306
    elif k == 6:
        return 407
    elif k == 5:
        return 58
    elif k == 4:
        return 11
    elif k == 3:
        return 3
    elif k == 2:
        return 1
    else:
        return None

def get_num_orbits_cum(k):
    if k < 2:
        return None
    elif k == 2:
        return get_num_orbits(k)
    else:
        return get_num_orbits(k) + get_num_orbits_cum(k - 1)

def calc_orbit_counts_autogen_graphlets(k):
    assert k in [6] # the orca method only works for k=6
    canon_list, orbit_map = read_in_canon_list_and_orbit_map(k)
    justk_orbit_counts = [None] * get_num_orbits(k)
    bvs = get_connected_bvs(canon_list)
    assert len(bvs) == get_num_graphlets(k)

    for bv in bvs:
        blantitl_el = get_bv_el_with_blantitl_orbit_nodes(bv, canon_list, orbit_map)
        blantitl_graph_path = get_tmp_path(f'graphlet_k{k}_bv{bv}.el')
        write_el_to_file(blantitl_el, blantitl_graph_path)
        p = run_orca_raw(5, blantitl_graph_path) # just run orca until 5, and run BLANT sample from 6 and up
        os.remove(blantitl_graph_path)
        orbit_lines = p.stdout.decode().strip().split('\n')[1:]

        for line in orbit_lines:
            splitted = line.split()
            node_name = splitted[0]
            blantitl_orbit_num = int(node_name[:-1])
            blantout_orbit_num = BLANTITL_TO_BLANTOUT_MAPPING[k][blantitl_orbit_num]
            orbits = splitted[1:]
            orbit_count = sum([1 if n != '0' else 0 for n in orbits])
            orbit_count += 1 # to include the orbit itself, since we'll only run blant or orca or whatever on k - 1

            if justk_orbit_counts[blantout_orbit_num] == None:
                justk_orbit_counts[blantout_orbit_num] = orbit_count
            else:
                assert justk_orbit_counts[blantout_orbit_num] == orbit_count

    # append justk to k5 one
    return CACHED_CUM_ORBIT_COUNTS[k - 1] + justk_orbit_counts
    
# directly use the graphlets directory which has the correct orbits
def calc_orbit_counts_direct(k):
    assert k in [4, 5]
    orbit_counts = [None] * get_num_orbits_cum(k)

    for graphlet_num in range(get_num_graphlets_cum(k)):
        p = run_orca_raw(k, get_base_graph_path(f'graphlets/graphlet{graphlet_num}'))
        orbit_lines = p.stdout.decode().strip().split('\n')[1:]

        for line in orbit_lines:
            splitted = line.split()
            node_name = splitted[0]
            orbit_num = int(node_name[:-1])
            orbits = splitted[1:]

            # we can't just sum all the orbits, we need to count how many are not zero (because if a node has a degree of 5 we only count that once for "an appearance of orbit 0)
            # we include the orbits effect on itself too
            orbit_count = sum([1 if n != '0' else 0 for n in orbits])

            if orbit_counts[orbit_num] == None:
                orbit_counts[orbit_num] = orbit_count
            else:
                assert orbit_counts[orbit_num] == orbit_count

        return orbit_counts
    
def calc_orbit_counts(k):
    USE_CACHE = True

    if USE_CACHE:
        if k in CACHED_CUM_ORBIT_COUNTS:
            return CACHED_CUM_ORBIT_COUNTS[k]
        else:
            return None
    else:
        if k in [6]:
            return calc_orbit_counts_autogen_graphlets(k)
        if k in [4, 5]:
            return calc_orbit_counts_direct(k)
        else:
            return None

def calc_weights(k):
    orbit_counts = calc_orbit_counts(k)
    weights = [1 - math.log(orbit_count) / math.log(get_num_orbits_cum(k)) for orbit_count in orbit_counts]
    return weights

def get_combined_odv_file(gtag, k, nstr, overwrite=True):
    cbodv_path = get_cbodv_path(gtag, k, nstr)

    if overwrite or not file_exists(cbodv_path):
        base_k = min(5, k - 1)
        assert base_k >= 4
        base_odv_dir = ODVDirectory(get_odv_path(gtag, base_k))
        assert k == base_k + 1 # since we hardcoded to read one splodv file for now
        spl_odv_dir = ODVDirectory(get_blantspl_path_nstr(gtag, k, nstr))
        assert base_odv_dir.get_nodes() == spl_odv_dir.get_nodes(), f'nodes not equal. base_odv_dir has {len(base_odv_dir.get_nodes())} nodes while spl_odv_dir has {len(spl_odv_dir.get_nodes())} nodes'
        out_odv = dict()

        for node in base_odv_dir.get_nodes():
            base_odv_list = base_odv_dir.get_odv(node).get_odv_list()
            spl_odv_list = spl_odv_dir.get_odv(node).get_odv_list()
            odv_list = base_odv_list + spl_odv_list
            assert len(odv_list) == get_num_orbits_cum(k)
            out_odv[node] = odv_list

        out_str = '\n'.join([f'{node} {" ".join(map(str, odv_list))}' for node, odv_list in out_odv.items()])
        write_to_file(out_str, cbodv_path)

    return ODVDirectory(cbodv_path)

class ODVDirectory:
    # file format: every line has node name, followed by orbit counts, separated by spaces
    # NODENAME 23 1 250 37 4 0 ...
    def __init__(self, fname):
        self._directory = dict()

        for line in open(fname, 'r'):
            line_split = line.strip().split()

            if len(line_split) == 1:
                continue # this is the initial first line of the .odv file that contains k
            
            node = line_split[0]
            odv_list = [int(s) for s in line_split[1:]]
            odv = ODV(node, odv_list)
            self._directory[node] = odv

    def get_odv(self, node):
        if node in self._directory:
            return self._directory[node]
        else:
            return None

    def get_nodes(self):
        return set(self._directory.keys())

    def __str__(self):
        return '\n'.join([f'{node}: {odv}' for node, odv in self._directory.items()])


class ODV:
    WEIGHTS = []
    WEIGHT_SUM = 0

    @staticmethod
    def set_weights_vars(k):
        ODV.WEIGHTS = calc_weights(k)
        ODV.WEIGHT_SUM = sum(ODV.WEIGHTS) # 45.08670802954777 <- calculated value from .sim file

    def __init__(self, node, odv_list):
        self._node = node
        self._odv_list = odv_list

    def get_similarity(self, other):
        if len(self._odv_list) == 0 or len(other._odv_list) == 0: # handle the case where the node is not connected to anything in one or both files, causing it to appear with no numbers after it in the .odv file
            return 0

        assert len(self._odv_list) == len(other._odv_list) == len(ODV.WEIGHTS), f'self: {len(self._odv_list)}, other: {len(other._odv_list)}, weights: {len(ODV.WEIGHTS)}, self._node: {self._node}, other._node: {other._node}'
        distance_sum = sum([self._get_single_orbit_similarity(m1, m2, i) for i, (m1, m2) in enumerate(zip(self._odv_list, other._odv_list))])
        weight_sum = ODV.WEIGHT_SUM
        return 1 - distance_sum / weight_sum

    def get_inequal_orbits(self, other):
        assert len(self._odv_list) == len(other._odv_list) == len(ODV.WEIGHTS), f'self: {len(self._odv_list)}, other: {len(other._odv_list)}, weights: {len(ODV.WEIGHTS)}'

        inequal_orbits = []

        for i, (o1, o2) in enumerate(zip(self._odv_list, other._odv_list)):
            if o1 != o2:
                inequal_orbits.append(i)

        return inequal_orbits

    def get_mean_similarity(self, other):
        return mean([self._get_single_orbit_mean_similarity(m1, m2) for m1, m2 in zip(self._odv_list, other._odv_list)])

    def get_odv_val(self, num):
        return self._odv_list[num]

    def get_odv_list(self):
        return self._odv_list

    def __str__(self):
        return ' '.join([str(n) for n in self._odv_list])

    @staticmethod
    def _get_single_orbit_mean_similarity(m1, m2):
        return 1 if m1 == m2 == 0 else min(m1, m2) / max(m1, m2)

    @staticmethod
    def _get_single_orbit_similarity(m1, m2, i):
        # the base of the log doesn't matter
        top_inner = math.log(m1 + 1) - math.log(m2 + 1)
        bot = math.log(max(m1, m2) + 2)
        return ODV.WEIGHTS[i] * abs(top_inner) / bot

def read_in_nodes_wo_deg1(gtag):
    graph_path = get_graph_path(gtag)
    nodes = read_in_nodes(graph_path)
    adj_set = read_in_adj_set(graph_path)
    new_nodes = [node for node in nodes if len(adj_set[node]) > 1]
    return new_nodes

def get_deg_sim(node1, node2, adj_set1, adj_set2, max_deg1, max_deg2):
    deg1 = len(adj_set1[node1])
    deg2 = len(adj_set2[node2])
    return (deg1 + deg2) / (max_deg1 + max_deg2)

# n is the number of orthologs to generate
# bn stands for "BLANT n" and is the number of samples
# if using bn, we're assuming that we're using a combined odv file over a normal odv file
def get_odv_orthologs(gtag1, gtag2, k, n, bnstr=None, no1=False, alpha=1):
    graph_path1 = get_graph_path(gtag1)
    graph_path2 = get_graph_path(gtag2)

    if no1:
        nodes1 = read_in_nodes_wo_deg1(gtag1)
        nodes2 = read_in_nodes_wo_deg1(gtag2)
    else:
        nodes1 = list(read_in_nodes(graph_path1))
        nodes2 = list(read_in_nodes(graph_path2))

    if bnstr == None:
        odv_path1 = get_odv_path(gtag1, k)
        odv_path2 = get_odv_path(gtag2, k)
    else:
        odv_path1 = get_cbodv_path(gtag1, k, bnstr)
        odv_path2 = get_cbodv_path(gtag2, k, bnstr)
        
    odv_dir1 = ODVDirectory(odv_path1)
    odv_dir2 = ODVDirectory(odv_path2)
    adj_set1 = read_in_adj_set(graph_path1)
    adj_set2 = read_in_adj_set(graph_path2)
    max_deg1 = get_max_deg(adj_set1)
    max_deg2 = get_max_deg(adj_set2)

    # assert n < len(nodes1) * len(nodes2), f'{n} must be >= {len(nodes1)} * {len(nodes2)}'

    top_n = [(-1, '', '')] * n
    heapq.heapify(top_n)
    tot_nodes = len(nodes1) # approximation for less incrementing
    proc_nodes = 0
    percent_printed = 0
    skip = 1

    # assert alpha == 1.0 # since we're commenting out deg_sim

    for node1 in nodes1:
        for i in range(0, len(nodes2), skip):
            node2 = nodes2[i]
            odv1 = odv_dir1.get_odv(node1)
            odv2 = odv_dir2.get_odv(node2)

            if odv1 == None or odv2 == None:
                continue

            odv_sim = odv1.get_similarity(odv2)
            deg_sim = get_deg_sim(node1, node2, adj_set1, adj_set2, max_deg1, max_deg2) # the reason we pass in max is so that we don't have to recalculate it every time we call this
            sim = alpha * odv_sim + (1 - alpha) * deg_sim

            # don't do min/max node just for sorting purposes, because the nodes come from two different graphs
            # min_node = min(node1, node2) BAD
            # max_node = max(node1, node2) BAD
            obj = (sim, node1, node2)
            min_top = heapq.heappushpop(top_n, obj)
        
        proc_nodes += 1

        if proc_nodes * 10000 / tot_nodes > percent_printed:
            percent_printed += 1
            print(f'{proc_nodes} / {tot_nodes}', file=sys.stderr)

    return sorted(top_n, reverse=True)

def analyze_mcl_test_data():
    nif1_path = get_data_path('mcl/mcl_test/ppi1.nif')
    nif2_path = get_data_path('mcl/mcl_test/ppi2.nif')
    ort_path = get_data_path('mcl/mcl_test/ppi1-ppi2.ort')
    ppi1_nodes = set()
    ppi2_nodes = set()

    with open(nif1_path, 'r') as nif1:
        for line in nif1:
            node1, node2, _ = line.strip().split('\t')
            ppi1_nodes.add(node1)
            ppi1_nodes.add(node2)

    with open(nif2_path, 'r') as nif2:
        for line in nif2:
            node1, node2, _ = line.strip().split('\t')
            ppi2_nodes.add(node1)
            ppi2_nodes.add(node2)

    ort_ppi1_nodes = set()
    ort_ppi2_nodes = set()

    with open(ort_path, 'r') as ort:
        for line in ort:
            print(line)
            node1, node2, _ = line.strip().split('\t')
            ort_ppi1_nodes.add(node1)
            ort_ppi2_nodes.add(node2)

    print(len(ppi1_nodes), len(ort_ppi1_nodes))
    print(len(ppi2_nodes), len(ort_ppi2_nodes))
    all_nodes = ppi1_nodes.union(ppi2_nodes)
    all_ort_nodes = ort_ppi1_nodes.union(ort_ppi2_nodes)
    print(len(all_nodes), len(all_ort_nodes))

def get_fake_ort_path(base, ext):
    return get_data_path(f'mcl/fake_ort/{base}.{ext}')

def get_odv_ort_path(gtag1, gtag2, k, n, bnstr=None, notes=''):
    base = f'{gtag1}-{gtag2}-k{k}-n{n}'

    if bnstr != None:
        base += f'-bn{bnstr}'
    
    if notes != '':
        base += f'-{notes}'

    return get_fake_ort_path(base, 'ort')

def get_default_odv_ort_path(gtag1, gtag2, notes=''):
    k = two_gtags_to_k(gtag1, gtag2)
    n = two_gtags_to_n(gtag1, gtag2)
    return get_odv_ort_path(gtag1, gtag2, k, n, notes=notes)

def read_in_odv_orts(path, include_score=True):
    from graph_helpers import unmark_node

    with open(path, 'r') as f:
        lines = f.readlines()
        splitted_strs = [line.strip().split('\t') for line in lines]

        if include_score:
            return [(unmark_node(node1), unmark_node(node2), float(score)) for node1, node2, score in splitted_strs]
        else:
            return [(unmark_node(node1), unmark_node(node2)) for node1, node2, _ in splitted_strs]

def odv_ort_file_to_nodes(path, left):
    from graph_helpers import unmark_node

    with open(path, 'r') as f:
        nodes = []

        for line in f:
            marked_node1, marked_node2, score = line.strip().split('\t')
            node1 = unmark_node(marked_node1)
            node2 = unmark_node(marked_node2)

            if left:
                nodes.append(node1)
            else:
                nodes.append(node2)

        return nodes

def odv_ort_to_nodes(odv_orts, left):
    nodes = list()

    for score, node1, node2 in odv_orts:
        if left:
            nodes.append(node1)
        else:
            nodes.append(node2)

    return nodes

def gen_fake_ort_from_sim(base, k, n):
    sim_path = get_fake_ort_path(base, 'sim')
    ort_path = get_fake_ort_path(f'{base}-{k}', 'ort')
    added_nodes = set()

    with open(sim_path, 'r') as sim_f:
        with open(ort_path, 'w') as ort_f:
            i = 0

            for line in sim_f:
                node1, node2, score = line.split()
                marked_node2 = f'sy05_{node2}'

                if i < n:
                    added_nodes.add(node1)
                    added_nodes.add(node2)
                    ort_f.write('\t'.join([node1, marked_node2, score]) + '\n')
                    i += 1
                else:
                    break

# function I used to validate the sim function based on Hayes' sim files
def validate_sim_function(gtag1, gtag2):
    FACTOR = 1_000_000

    odv_path1 = get_odv_path(gtag1, 5)
    odv_path2 = get_odv_path(gtag2, 5)
    odv_dir1 = ODVDirectory(odv_path1)
    odv_dir2 = ODVDirectory(odv_path2)
    sim_path = get_fake_ort_path(f'{gtag1}-{gtag2}', 'sim')

    tot_diff = 0
    tot_pairs = 0
    num_gt10 = 0

    with open(sim_path, 'r') as sim_file:
        for line in sim_file:
            node1, node2, sim = line.strip().split()
            sim = float(sim)
            sim_non_decimal = int(sim * FACTOR) # cuz the sim_path values are rounded to six
            odv1 = odv_dir1.get_odv(node1)
            odv2 = odv_dir2.get_odv(node2)
            my_sim = odv1.get_similarity(odv2)
            my_sim_non_decimal = int(my_sim * FACTOR)
            tot_diff += abs(sim_non_decimal - my_sim_non_decimal)
            tot_pairs += 1

            if tot_pairs % 5000 == 0:
                print(tot_pairs, '/', 1004 ** 2)

            if tot_pairs > 10000:
                break

    avg_diff = (tot_diff / tot_pairs) / FACTOR
    print(f'avg_diff: {avg_diff}')
    print(f'num_gt10: {num_gt10}')

def odv_ort_to_str(odv_ort, mark1, mark2):
    return '\n'.join([f'{mark1}_{node1}\t{mark2}_{node2}\t{score}' for score, node1, node2 in odv_ort])

def make_odv_ort_1to1(odv_ort):
    odv_ort_1to1 = []
    used_nodes1 = set()
    used_nodes2 = set()

    for node1, node2 in odv_ort:
        if node1 in used_nodes1 or node2 in used_nodes2:
            continue

        odv_ort_1to1.append((node1, node2))
        used_nodes1.add(node1)
        used_nodes2.add(node2)

    return odv_ort_1to1

def get_odv_alignment(odv_ort, adj_set1, adj_set2):
    from analysis_helpers import get_s3
    
    alignment = []
    step_size = len(odv_ort) // 10
    n = step_size

    while n <= len(odv_ort):
        alignment = odv_ort[:n]
        s3 = get_s3(alignment, adj_set1, adj_set2)

        if s3 < 0.8:
            break
        
        n += step_size

    return alignment


if __name__ == '__main__':
    from graph_helpers import get_graph_path, read_in_adj_set
    from analysis_helpers import get_deg_distr, print_deg_distr

    gtag1 = sys.argv[1]
    gtag2 = sys.argv[2]
    k = two_gtags_to_k(gtag1, gtag2)
    ODV.set_weights_vars(k)
    n = two_gtags_to_n(gtag1, gtag2)
    odv_ort = get_odv_orthologs(gtag1, gtag2, k, n)
    print(odv_ort_to_str(odv_ort, '', ''))
