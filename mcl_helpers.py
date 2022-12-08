#!/pkg/python/3.7.4/bin/python3
import sys
import os
from collections import defaultdict

def read_in_slashes_m2m(m2m_path):
    from graph_helpers import unmark_node
    
    m2m_pairs = []

    with open(m2m_path, 'r') as f:
        for alignment in f:
            for line in alignment.strip().split('\t'):
                node1, node2 = line.split('/')
                m2m_pairs.append((unmark_node(node1), unmark_node(node2)))

    return m2m_pairs

def read_in_slashes_alignments(path):
    from graph_helpers import unmark_node
    
    alignments = []

    with open(path, 'r') as f:
        for alignment in f:
            new_alignment = []

            for line in alignment.strip().split('\t'):
                node1, node2 = line.split('/')
                new_alignment.append((unmark_node(node1), unmark_node(node2)))

            alignments.append(new_alignment)

    return alignments

def get_nif_str(el):
    return('\n'.join([f'{node1}\t{node2}\t1' for node1, node2 in el]))

def gen_nif_file(gtag, overwrite=False):
    from graph_helpers import get_nif_path, get_marked_el
    from file_helpers import file_exists, write_to_file

    nif_path = get_nif_path(gtag)

    if overwrite or not file_exists(nif_path):
        marked_el = get_marked_el(gtag)
        write_to_file(get_nif_str(marked_el), nif_path)
    else:
        print(f'using old nif file for {gtag}', file=sys.stderr)

def gen_odv_ort_file(gtag1, gtag2, override_k=None, bnstr=None, overwrite=False, notes=''):
    from graph_helpers import gtag_to_mark, read_in_adj_set, get_graph_path
    from odv_helpers import get_odv_orthologs, odv_ort_to_str, get_odv_ort_path, two_gtags_to_k, two_gtags_to_n, ODV, odv_ort_to_nodes
    from file_helpers import file_exists, write_to_file
    from analysis_helpers import get_deg_distr

    k = two_gtags_to_k(gtag1, gtag2, override_k=override_k)
    n = two_gtags_to_n(gtag1, gtag2)
    ODV.set_weights_vars(k)
    ort_path = get_odv_ort_path(gtag1, gtag2, k, n, bnstr=bnstr, notes=notes)

    if overwrite or not file_exists(ort_path):
        if notes == 'no1':
            no1 = True
            alpha = 1
        elif notes == 'a80':
            no1 = False
            alpha = 0.8
        elif notes == 'norm':
            no1 = False
            alpha = 1
            
        odv_ort = get_odv_orthologs(gtag1, gtag2, k, n, bnstr=bnstr, no1=no1, alpha=alpha)
        mark1 = gtag_to_mark(gtag1)
        mark2 = gtag_to_mark(gtag2)
        ort_str = odv_ort_to_str(odv_ort, mark1, mark2)
        write_to_file(ort_str, ort_path)
    else:
        print(f'using old odv ort file for {gtag1}-{gtag2}', file=sys.stderr)

def copy_to_out(gtag1, gtag2, notes=''):
    from bash_helpers import run_outsend
    from graph_helpers import get_nif_path
    from odv_helpers import get_odv_ort_path, two_gtags_to_k, two_gtags_to_n

    nif_path1 = get_nif_path(gtag1)
    nif_path2 = get_nif_path(gtag2)
    k = two_gtags_to_k(gtag1, gtag2)
    n = two_gtags_to_n(gtag1, gtag2)
    ort_path = get_odv_ort_path(gtag1, gtag2, k, n, notes=notes)
    run_outsend(nif_path1)
    run_outsend(nif_path2)
    run_outsend(ort_path)

def get_mcl_out_fname(gtag1, gtag2, k, n, notes=''):
    base = f'{gtag1}-{gtag2}-k{k}-n{n}'

    if notes != '':
        base += f'-{notes}'

    return f'{base}.txt'

def read_in_mcl_out(path):
    with open(path, 'r') as f:
        pass

def get_mcl_out_path(gtag1, gtag2, k, n, notes=''):
    from file_helpers import get_data_path
    return get_data_path(f'mcl/{get_mcl_out_fname(gtag1, gtag2, k, n, notes=notes)}')

def get_mcl_paths(gtag1, gtag2, k, n, notes=''):
    out_path = get_mcl_out_path(gtag1, gtag2, k, n, notes=notes)
    out_path_base = '.'.join(out_path.split('.')[:-1])
    ag_path = out_path_base + '.ag'
    time_path = out_path_base + '.time'

    return out_path, ag_path, time_path

def take_from_out(gtag1, gtag2, notes=''):
    from bash_helpers import run_outtake
    from odv_helpers import two_gtags_to_k, two_gtags_to_n

    k = two_gtags_to_k(gtag1, gtag2)
    n = two_gtags_to_n(gtag1, gtag2)
    out_fname = get_mcl_out_fname(gtag1, gtag2, k, n, notes=notes)
    here_path = get_mcl_out_path(gtag1, gtag2, k, n, notes=notes)
    run_outtake(out_fname, here_path) # note: the out here means the output file of mcl, while the out param of run_outtake refers to the remote git repo called "out" (which is named vmcopy lol)

def process_mcl(gtag1, gtag2, notes=''):
    from odv_helpers import two_gtags_to_k, two_gtags_to_n
    from analysis_helpers import print_distr, get_s3, get_alignment_nc
    from graph_helpers import get_graph_path, read_in_adj_set
    from ortholog_helpers import get_g1_to_g2_orthologs

    k = two_gtags_to_k(gtag1, gtag2)
    n = two_gtags_to_n(gtag1, gtag2)
    adj_set1 = read_in_adj_set(get_graph_path(gtag1))
    adj_set2 = read_in_adj_set(get_graph_path(gtag2))
    g1_to_g2_ort = get_g1_to_g2_orthologs(gtag1, gtag2)
    out_path = get_mcl_out_path(gtag1, gtag2, k, n, notes=notes)
    print(f'evaluating {out_path}')
    alignments = read_in_slashes_alignments(out_path)

    if False:
        lens = defaultdict(int)
            
        for alignment in alignments:
            lens[len(alignment)] += 1
    
        print_distr(lens, 'alignment size')

    if False:
        for alignment in alignments:
            alignment_size = len(alignment)
            s3 = get_s3(alignment, adj_set1, adj_set2)

            if alignment_size == 1 or s3 == 0:
                continue
            
            print(alignment_size, s3)

    if False:
        for alignment in alignments:
            alignment_size = len(alignment)
            s3 = get_s3(alignment, adj_set1, adj_set2)
            nc = get_alignment_nc(alignment, g1_to_g2_ort)

            if alignment_size < 8 or s3 == 0:
                continue
            
            print(s3, nc)

def get_mcl_tfp_stats(alignments, g1_to_g2_ort, adj_set1, adj_set2):
    from analysis_helpers import get_topofunc_perfect_alignments
    
    alignments = [alignment for alignment in alignments if len(alignment) >= 8]
    tfp_aligns = get_topofunc_perfect_alignments(alignments, g1_to_g2_ort, adj_set1, adj_set2)
    lens = [len(align) for align in tfp_aligns]
    max_len = max(lens) if len(lens) > 0 else 0
    return len(tfp_aligns), max_len
        
def prepare_mcl(gtag1, gtag2, notes=''):
    from bash_helpers import run_orca_for_gtag

    gen_nif_file(gtag1)
    gen_nif_file(gtag2)
    run_orca_for_gtag(gtag1)
    run_orca_for_gtag(gtag2)
    gen_odv_ort_file(gtag1, gtag2, notes=notes)

def full_local_run_mcl(gtag1, gtag2, notes=''):
    from bash_helpers import run_align_mcl

    prepare_mcl(gtag1, gtag2, notes=notes)
    p = run_align_mcl(gtag1, gtag2, notes=notes)

    if p != None:
        p.wait()

    process_mcl(gtag1, gtag2, notes=notes)

def wayne_copy_mcl(gtag1, gtag2, notes=''):
    from bash_helpers import run_cmd
    from general_helpers import get_wayne_path
    from odv_helpers import two_gtags_to_k, two_gtags_to_n

    k = two_gtags_to_k(gtag1, gtag2)
    n = two_gtags_to_n(gtag1, gtag2)
    old_path = get_mcl_out_path(gtag1, gtag2, k, n, notes=notes)
    new_path = get_wayne_path(f'mcl/{gtag1}-{gtag2}-mclseeds.txt')
    p = run_cmd(f'cp {old_path} {new_path}')

def clean_mcl_single(gtag, notes=''):
    from odv_helpers import gtag_to_k, get_odv_path
    from graph_helpers import get_nif_path

    k = gtag_to_k(gtag)
    nif_path = get_nif_path(gtag)
    odv_path = get_odv_path(gtag, k)

    for path in [nif_path, odv_path]:
        try:
            os.remove(path)
            print(f'removed {path}')
        except Exception as e:
            print(f'failed to remove {path} because of {e}')
            
def clean_mcl_pair(gtag1, gtag2, notes=''):
    from odv_helpers import two_gtags_to_k, two_gtags_to_n, get_odv_ort_path
    
    k = two_gtags_to_k(gtag1, gtag2)
    n = two_gtags_to_n(gtag1, gtag2)
    ort_path = get_odv_ort_path(gtag1, gtag2, k, n, notes=notes)
    out_path, ag_path, time_path = get_mcl_paths(gtag1, gtag2, k, n, notes=notes)

    for path in [ort_path, mcl_out_path, ag_path, time_path]:
        try:
            os.remove(path)
            print(f'removed {path}')
        except Exception as e:
            print(f'failed to remove {path} because of {e}')

def clean_mcl(gtag1, gtag2, notes=''):
    clean_mcl_single(gtag1, notes=notes)
    clean_mcl_single(gtag2, notes=notes)
    clean_mcl_pair(gtag1, gtag2, notes=notes)
    
if __name__ == '__main__':
    mode = sys.argv[1]
    gtag1 = sys.argv[2]
    gtag2 = sys.argv[3]
    notes = sys.argv[4] if len(sys.argv) > 4 else ''

    if mode == 'prep':
        prepare_mcl(gtag1, gtag2, notes=notes)
        # copy_to_out(gtag1, gtag2, notes=notes)
    elif mode == 'proc':
        # take_from_out(gtag1, gtag2, notes=notes)
        process_mcl(gtag1, gtag2, notes=notes)
    elif mode == 'full':
        full_local_run_mcl(gtag1, gtag2, notes=notes)
    elif mode == 'clean':
        clean_mcl(gtag1, gtag2, notes=notes)
    elif mode == 'wcpy':
        wayne_copy_mcl(gtag1, gtag2, notes=notes)
    else:
        print('USAGE: mcl_helpers.py gtag1 gtag2 mode, where mode is prep or proc')
