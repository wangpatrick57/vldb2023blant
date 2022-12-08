#!/pkg/python/3.7.4/bin/python3
from graph_helpers import *
from index_helpers import *
from file_helpers import *
import subprocess
import time
import sys

def bool_conv(b):
    return 1 if b else 0

def run_blant(gtag, lDEG=2, alph=True, algo='stairs', overwrite=False):
    assert alph != None # alph can't be None because we need a different .sh script for that
    graph_path = get_graph_path(gtag)
    out_path = get_index_path(gtag, lDEG=lDEG, alph=alph, algo=algo)

    if overwrite or not file_exists(out_path):
        common_settings = f'{graph_path} {lDEG} {bool_conv(alph)} {out_path}'

        if algo == None:
            cmd = f'run_blant_default_custom.sh {common_settings}'
        else:
            cmd = f'run_blant_default_custom_algo.sh {algo} {common_settings}'

        p = subprocess.Popen(cmd.split())
    else:
        p = None
        print(f'using old index file for {gtag}', file=sys.stderr)

    return p, out_path

def run_blant_sample_raw(k, n, graph_path):
    out_path = get_tmp_path(f'blantspl-k{k}-n{n}.out')
    # we need to call run_blant_sample.sh because it takes care of cd for us
    cmd = f'run_blant_sample.sh {graph_path} {k} {n} {out_path}'
    p = subprocess.run(cmd.split(), capture_output=True)

    with open(out_path, 'r') as f:
        out = ''.join(f.readlines())

    os.remove(out_path)
    return out

# should be modified to call run_blant_sample_raw to be consistent with orca, but I'm too lazy to do so
def run_blant_sample(gtag, k, n, overwrite=False):
    from odv_helpers import get_blantspl_path
    
    graph_path = get_graph_path(gtag)
    out_path = get_blantspl_path(gtag, k, n)

    if overwrite or not file_exists(out_path):
        cmd = f'run_blant_sample.sh {graph_path} {k} {n} {out_path}'
        p = subprocess.Popen(cmd.split())
    else:
        p = None
        print(f'using old index file for {gtag}', file=sys.stderr)

    return p, out_path

def run_outsend(path):
    cmd = f'outsend {path}'
    subprocess.run(cmd.split())

def run_outtake(out_fname, here_path):
    cmd = f'outtake {out_fname} {here_path}'
    subprocess.run(cmd.split())

def run_cmd(cmd):
    p = subprocess.run(cmd.split())
    return p

def run_orca_raw(k, graph_path):
    cmd = f'orca.sh {k} {graph_path}'
    p = subprocess.run(cmd.split(), capture_output=True)
    return p

def run_orca_for_gtag(gtag, override_k=None, overwrite=False):
    from graph_helpers import get_graph_path
    from file_helpers import write_to_file, file_exists
    from odv_helpers import get_odv_path, gtag_to_k
    
    k = gtag_to_k(gtag, override_k=override_k)
    odv_path = get_odv_path(gtag, k)

    if overwrite or not file_exists(odv_path):
        graph_path = get_graph_path(gtag)
        p = run_orca_raw(k, graph_path)
        out = p.stdout.decode('utf-8')
        write_to_file(out, odv_path)
    else:
        print(f'using old odv file for {gtag}', file=sys.stderr)

def run_align_mcl(gtag1, gtag2, overwrite=False, notes=''):
    from graph_helpers import get_nif_path
    from odv_helpers import two_gtags_to_k, two_gtags_to_n, get_odv_ort_path
    from mcl_helpers import get_mcl_out_path
    from file_helpers import file_exists, get_home_path, remove_extension
    MCL_SCRIPT_PATH = get_home_path('alignMCL/raw_run_align_mcl.sh')

    k = two_gtags_to_k(gtag1, gtag2)
    n = two_gtags_to_n(gtag1, gtag2)
    out_path = get_mcl_out_path(gtag1, gtag2, k, n, notes=f'{notes}')

    if overwrite or not file_exists(out_path):
        out_arg = remove_extension(out_path)
        nif_path1 = get_nif_path(gtag1)
        nif_path2 = get_nif_path(gtag2)
        odv_ort_path = get_odv_ort_path(gtag1, gtag2, k, n, notes=notes)

        from file_helpers import get_num_lines
        print(out_path, nif_path1, nif_path2, odv_ort_path, sep='\n')
        cmd = f'{MCL_SCRIPT_PATH} {nif_path1} {nif_path2} {odv_ort_path} {out_arg}'
        p = subprocess.Popen(cmd.split())
    else:
        p = None
        print(f'using old mcl out file for {gtag1}-{gtag2}', file=sys.stderr)

    return p

if __name__ == '__main__':
    from odv_helpers import two_gtags_to_n

    for gtag in get_tprl_gtags():
        print(f'starting {gtag}')
        run_orca_for_gtag(gtag)
        print(f'done with {gtag}')
