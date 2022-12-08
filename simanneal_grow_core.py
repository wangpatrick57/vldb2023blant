#!/bin/python3
# start with random seed
# at every step, pick a random seed to flip the bit of (remove if it's in the alignment, add if it's not)
# energy function: number of nodes * s3 ** 2
# disallow moves that send the s3 below the threshold
from all_helpers import *
import sys
import copy
import unittest
import random
import curses

STANDARD_P_FUNC = 'standard_p_func'
ALWAYS_P_FUNC = 'always_p_func'

# blocks are the building block alignments
def seeds_to_blocks(seeds):
    seeds.sort()
    blocks = []

    for sid, nodes1, nodes2 in seeds:
        blocks.append(list(zip(nodes1, nodes2)))

    return blocks

def get_mcl_blocks(gtag1, gtag2, min_size, min_s3, notes='no1'):
    adj_set1 = read_in_adj_set(get_graph_path(gtag1))
    adj_set2 = read_in_adj_set(get_graph_path(gtag2))
    out_path = get_mcl_out_path(gtag1, gtag2, two_gtags_to_k(gtag1, gtag2), two_gtags_to_n(gtag1, gtag2), notes='no1')
    alignments = read_in_slashes_alignments(out_path)
    alignments = get_clean_alignments(alignments, adj_set1, adj_set2)
    return [alignment for alignment in alignments if len(alignment) >= min_size and get_s3(alignment, adj_set1, adj_set2) >= min_s3]

def energy(size):
    return -size

# bad means increase so it has to be positive
def get_temperature_endpoints(avg_bad_energy_delta):
    assert avg_bad_energy_delta >= 1 # just to keep the scale sane

    start_t = 1

    while math.exp(-avg_bad_energy_delta / start_t) < 0.99:
        start_t *= 10

    while math.exp(-avg_bad_energy_delta / start_t) > 0.99:
        start_t /= 2

    end_t = 1

    while math.exp(-avg_bad_energy_delta / end_t) > 0.000001:
        end_t /= 10

    while math.exp(-avg_bad_energy_delta / end_t) < 0.000001:
        end_t *= 2

    return start_t, end_t
        
def s3_frac_to_s3(s3_frac):
    if s3_frac[1] == 0:
        return None
    else:
        return s3_frac[0] / s3_frac[1]

class SimAnnealGrow:
    # blocks are the building block alignments
    # WHEN TO PUT PARAMS IN INIT AND WHEN TO PUT THEM IN RUN
    # if it's used only in run, put it in run
    # if it's used in functions other than run, put it in init
    # write_progress_info really should go in run but that's a todo for later
    def __init__(self, blocks, adj_set1, adj_set2, p_func=STANDARD_P_FUNC, s3_threshold=1, write_progress_info=None, read_use_block_path=None):
        assert all(is_well_formed_alignment(block) for block in blocks)
        assert is_symmetric_adj_set(adj_set1)
        assert is_symmetric_adj_set(adj_set2)
        self._blocks = blocks
        self._adj_set1 = adj_set1
        self._adj_set2 = adj_set2
        self._p_func = p_func
        self._s3_threshold = s3_threshold
        self._read_use_block_path = read_use_block_path

        if write_progress_info != None:
            assert type(write_progress_info) is tuple and len(write_progress_info) == 2
            assert type(write_progress_info[0]) is str
            assert type(write_progress_info[1]) is int
            self._write_progress = True
            self._progress_path = write_progress_info[0]
            self._write_every_k = write_progress_info[1]

            with open(self._progress_path, 'w') as f:
                pass
        else:
            self._write_progress = False
        
        self._pair_multiset = dict()
        self._forward_mapping = dict()
        self._reverse_mapping = dict()
        self._s3_frac = [0, 0]
        # store numer and denom for s3
        self._reset()

    # returns a valid move (without caring about the energy change)
    def _neighbor(self):
        block_i = self._get_random_block_i()
        
        while not self._move_is_valid(block_i):
            block_i = self._get_random_block_i()

        return block_i

    def _p(self, e, e_new, temperature):
        if self._p_func == STANDARD_P_FUNC:
            if e_new < e:
                return 1.0
            else:
                return 0 if temperature == 0 else math.exp(-(e_new - e) / temperature)
        elif self._p_func == ALWAYS_P_FUNC:
            return 1.0
        else:
            raise AssertionError()
    
    def _move_is_valid(self, block_i):
        is_adding = not self._use_block[block_i]
        block = self._blocks[block_i]
        
        # we can't make the mapping non-injective with an add
        if is_adding:
            for node1, node2 in block:
                if node1 in self._forward_mapping:
                    if self._forward_mapping[node1] != node2:
                        return False

                if node2 in self._reverse_mapping:
                    if self._reverse_mapping[node2] != node1:
                        return False

        # we can't go below the threshold
        updated_s3_frac = self._get_updated_s3_frac(block_i)
        updated_s3 = s3_frac_to_s3(updated_s3_frac)
        
        # even the first move has to have good enough s3, so we don't have a special case for that. moves add blocks at a time not nodes, so we don't treat the first move differently
        if self._s3_threshold == None:
            pass # if the threshold is None, we'll even allow moves that give undefined s3 scores
        else:
            is_removing_last_block = not is_adding and self._num_used_blocks == 1

            if not is_removing_last_block:
                if updated_s3 == None or updated_s3 < self._s3_threshold:
                    return False

        return True

    def _get_random_block_i(self):
        return random.randrange(len(self._blocks))
    
    # call this to modify use block in order to do some bookkeeping
    def _make_move(self, block_i):
        assert self._move_is_valid(block_i)

        # s3 needs to be updated before the bookkeeping stuff
        self._s3_frac = self._get_updated_s3_frac(block_i)

        # update bookkeeping stuff
        is_adding = not self._use_block[block_i]
        
        if is_adding:
            self._update_block_bookkeeping_after_add(block_i)
        else:
            self._update_block_bookkeeping_after_remove(block_i)

        # modify use block last since everything else depends on this value being old to calculate is_adding correctly
        self._use_block[block_i] = not self._use_block[block_i]

    def _update_block_bookkeeping_after_add(self, block_i):
        block = self._blocks[block_i]
        self._num_used_blocks += 1
                    
        for node1, node2 in block:
            pair = (node1, node2)

            if pair not in self._pair_multiset:
                self._pair_multiset[pair] = 1
                self._forward_mapping[node1] = node2
                self._reverse_mapping[node2] = node1
            else:
                self._pair_multiset[pair] += 1

    def _update_block_bookkeeping_after_remove(self, block_i):
        block = self._blocks[block_i]
        self._num_used_blocks -= 1

        for node1, node2 in block:
            pair = (node1, node2)
            self._pair_multiset[pair] -= 1

            if self._pair_multiset[pair] == 0:
                del self._pair_multiset[pair]
                del self._forward_mapping[node1]
                del self._reverse_mapping[node2]
                
    def _get_updated_s3_frac(self, block_i):
        # this function does not update any internal fields
        block = self._blocks[block_i]
        is_adding = not self._use_block[block_i]
        delta_pairs = self._get_delta_pairs(block_i)
        curr_aligned_pairs = set(self._pair_multiset.keys())
        delta = 1 if is_adding else -1
        updated_s3_frac = copy.copy(self._s3_frac)

        for node1, node2 in delta_pairs:
            pair = (node1, node2)
            
            for curr_node1, curr_node2 in curr_aligned_pairs:
                has_edge1 = node1 in self._adj_set1[curr_node1]
                has_edge2 = node2 in self._adj_set2[curr_node2]

                if has_edge1 and has_edge2:
                    updated_s3_frac[0] += delta
                
                if has_edge1 or has_edge2:
                    updated_s3_frac[1] += delta
                    
            if is_adding:
                curr_aligned_pairs.add(pair)
            else:
                curr_aligned_pairs.remove(pair)

        return updated_s3_frac
    
    def _reset(self):
        self._use_block = [False] * len(self._blocks)
        self._num_used_blocks = 0
        
        if self._read_use_block_path != None:
            read_use_block = []
            
            with open(self._read_use_block_path) as f:
                for line in f:
                    read_use_block.append(bool(int(line.strip())))

            assert len(read_use_block) == len(self._use_block)

            for block_i, use_block in enumerate(read_use_block):
                if use_block:
                    print(f'attempting {block_i}')
                    self._make_move(block_i)
                    print(f'did {block_i}')
            
    def _get_delta_pairs(self, block_i):
        is_adding = not self._use_block[block_i]
        block = self._blocks[block_i]

        if is_adding:
            return {pair for pair in block if pair not in self._pair_multiset}
        else:
            return {pair for pair in block if self._pair_multiset[pair] == 1}
        
    def _get_num_delta_pairs(self, block_i):
        is_adding = not self._use_block[block_i]
        abs_num_delta_pairs = len(self._get_delta_pairs(block_i))

        if is_adding:
            return abs_num_delta_pairs
        else:
            return -1 * abs_num_delta_pairs
        
    def run(self, k_max=None, auto_k=None, silent=False, write_use_block_path=None):
        if k_max != None:
            assert auto_k == None
            assert type(k_max) is int
        elif auto_k != None:
            assert k_max == None
            assert type(auto_k) is tuple
            assert len(auto_k) == 2
            assert type(auto_k[0]) == int
            assert type(auto_k[1]) == float
            auto_k_window = auto_k[0]
            auto_k_inc_ratio = auto_k[1]
            auto_k_granularity = 100
            assert auto_k[0] % auto_k_granularity == 0
            back_shift = auto_k_window // auto_k_granularity
            auto_k_size_log = []
        else:
            raise AssertionError('either k_max or auto_k has to be set')

        self._reset()
        
        lowest_energy = None
        size_after_last_move = 0
        best_alignment = []
        best_use_block = self._use_block
        start_t, end_t = get_temperature_endpoints(10)
        start_time = time.time()
        k = 0
        inc_ratio = None

        while True:
            if not silent:
                if k % 1000 == 0:
                    if k_max != None:
                        print(f'{k} / {k_max}, {time.time() - start_time:.2f}s', file=sys.stderr)
                    else:
                        print(f'{k} iter, {inc_ratio} ratio, {time.time() - start_time:.2f}s', file=sys.stderr)

            if k_max != None:
                if k >= k_max:
                    break
            elif auto_k != None:
                if k % auto_k_granularity == 0:
                    auto_k_size_log.append(len(best_alignment))
                    curr_log_i = k // auto_k_granularity

                    if len(auto_k_size_log) > back_shift:
                        curr_size = auto_k_size_log[curr_log_i]
                        old_size = auto_k_size_log[curr_log_i - back_shift]
                        inc_ratio = None if old_size == 0 else curr_size / old_size - 1

                        if inc_ratio != None and inc_ratio < auto_k_inc_ratio:
                            break

            if k_max != None:
                temperature = (1 - k / (k_max - 1)) * (start_t - end_t) + end_t
            else:
                temperature = 1
                
            next_block_i = self._neighbor()
            num_pairs = len(self._pair_multiset)
            delta_pairs = self._get_num_delta_pairs(next_block_i)
            curr_energy = energy(num_pairs)
            new_energy = energy(num_pairs + delta_pairs)

            if self._p(curr_energy, new_energy, temperature) >= random.random():
                self._make_move(next_block_i)
                size_after_last_move = len(self._pair_multiset)

                if lowest_energy == None or new_energy < lowest_energy:
                    lowest_energy = new_energy
                    best_alignment = list(self._pair_multiset)
                    best_use_block = list(self._use_block)

            # write regardless of whether we made a move or not
            if self._write_progress and k % self._write_every_k == 0:
                with open(self._progress_path, 'a') as f:
                    f.write(f'{k}\t{size_after_last_move}\t{len(best_alignment)}\n')
            k += 1

        if write_use_block_path != None:
            write_to_file('\n'.join(['1' if ub else '0' for ub in self._use_block]), write_use_block_path)
            
        return best_alignment

class TestSimAnnealGrow(unittest.TestCase):
    def setUp(self):
        self.blocks_sagrow = self.get_blocks_sagrow()
        self.s3_sagrow = self.get_s3_sagrow()
        self.blank_sagrow = self.get_blank_sagrow()
        
    def test_null_params(self):
        try:
            sagrow = SimAnnealGrow([], dict(), dict())
            found_error = False
        except:
            found_error = True

        self.assertFalse(found_error, 'Null params didn\'t get through')
        
    def test_well_formed_alignment_assertion(self):
        try:
            sagrow = SimAnnealGrow([('a', 'b'), ('a', 'b', 'c')], dict(), dict())
            found_error = False
        except:
            found_error = True

        self.assertTrue(found_error, 'Non-well-formed alignment got through')

    def test_symmetric_adj_set_assertion(self):
        try:
            sagrow = SimAnnealGrow([], {'a': {'b'}}, dict())
            found_error = False
        except:
            found_error = True

        self.assertTrue(found_error, 'Asymmetric adj_set got through')

        try:
            sagrow = SimAnnealGrow([], dict(), {'a': {'b', 'c'}, 'b': {'c'}, 'c': {'a'}})
            found_error = False
        except:
            found_error = True

        self.assertTrue(found_error, 'Asymmetric adj_set got through')
        
    def get_blocks_sagrow(self):
        # used to test block behavior
        return SimAnnealGrow([
            [('a', 'B')],
            [('b', 'A')],
            [('c', 'D'), ('a', 'B')],
            [('c', 'D'), ('e', 'F')],
            [('a', 'X')],
            [('y', 'B')]
        ], {
            'a': {},
            'b': {},
            'c': {},
            'e': {},
            'y': {},
        }, {
            'B': {},
            'A': {},
            'D': {},
            'F': {},
            'X': {},
        }, s3_threshold=None) # set threshold to None so all moves are valid

    def test_non_injective_add(self):
        self.assertTrue(self.blocks_sagrow._move_is_valid(4))
        self.assertTrue(self.blocks_sagrow._move_is_valid(5))
        self.blocks_sagrow._make_move(0)
        self.assertFalse(self.blocks_sagrow._move_is_valid(4))
        self.assertFalse(self.blocks_sagrow._move_is_valid(5))

    def test_add_get_num_delta_pairs(self):
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(0), 1)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(1), 1)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(2), 2)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(3), 2)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(4), 1)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(5), 1)
        self.blocks_sagrow._make_move(0)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(1), 1)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(2), 1)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(3), 2)
        self.blocks_sagrow._make_move(3)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(1), 1)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(2), 0)

    def test_remove_get_num_delta_pairs(self):
        self.blocks_sagrow._make_move(0)
        self.blocks_sagrow._make_move(1)
        self.blocks_sagrow._make_move(2)
        self.blocks_sagrow._make_move(3)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(0), 0)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(1), -1)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(2), 0)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(3), -1)
        self.blocks_sagrow._make_move(1)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(0), 0)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(2), 0)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(3), -1)
        self.blocks_sagrow._make_move(0)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(2), -1)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(3), -1)
        self.blocks_sagrow._make_move(3)
        self.assertEqual(self.blocks_sagrow._get_num_delta_pairs(2), -2)
        
    def test_add_bookkeeping_update(self):
        self.blocks_sagrow._make_move(0)
        self.assertEqual(self.blocks_sagrow._pair_multiset, {('a', 'B'): 1})
        self.assertEqual(self.blocks_sagrow._forward_mapping, {'a': 'B'})
        self.assertEqual(self.blocks_sagrow._reverse_mapping, {'B': 'a'})
        self.blocks_sagrow._make_move(2)
        self.assertEqual(self.blocks_sagrow._pair_multiset, {('a', 'B'): 2, ('c', 'D'): 1})
        self.assertEqual(self.blocks_sagrow._forward_mapping, {'a': 'B', 'c': 'D'})
        self.assertEqual(self.blocks_sagrow._reverse_mapping, {'B': 'a', 'D': 'c'})
        self.blocks_sagrow._make_move(3)
        self.assertEqual(self.blocks_sagrow._pair_multiset, {('a', 'B'): 2, ('c', 'D'): 2, ('e', 'F'): 1})
        self.assertEqual(self.blocks_sagrow._forward_mapping, {'a': 'B', 'c': 'D', 'e': 'F'})
        self.assertEqual(self.blocks_sagrow._reverse_mapping, {'B': 'a', 'D': 'c', 'F': 'e'})

    def test_remove_bookkeeping_update(self):
        self.blocks_sagrow._make_move(0)
        self.blocks_sagrow._make_move(2)
        self.blocks_sagrow._make_move(3)
        self.assertEqual(self.blocks_sagrow._pair_multiset, {('a', 'B'): 2, ('c', 'D'): 2, ('e', 'F'): 1})
        self.assertEqual(self.blocks_sagrow._forward_mapping, {'a': 'B', 'c': 'D', 'e': 'F'})
        self.assertEqual(self.blocks_sagrow._reverse_mapping, {'B': 'a', 'D': 'c', 'F': 'e'})
        self.blocks_sagrow._make_move(2)
        self.assertEqual(self.blocks_sagrow._pair_multiset, {('a', 'B'): 1, ('c', 'D'): 1, ('e', 'F'): 1})
        self.assertEqual(self.blocks_sagrow._forward_mapping, {'a': 'B', 'c': 'D', 'e': 'F'})
        self.assertEqual(self.blocks_sagrow._reverse_mapping, {'B': 'a', 'D': 'c', 'F': 'e'})
        self.blocks_sagrow._make_move(3)
        self.assertEqual(self.blocks_sagrow._pair_multiset, {('a', 'B'): 1})
        self.assertEqual(self.blocks_sagrow._forward_mapping, {'a': 'B'})
        self.assertEqual(self.blocks_sagrow._reverse_mapping, {'B': 'a'})

    def get_s3_sagrow(self):
        return SimAnnealGrow([
            [('a', 'A'), ('b', 'B')],
            [('c', 'C')],
            [('d', 'D')],
            [('e', 'E')],
            [('a', 'A'), ('c', 'C')],
        ], {
            'a': {'b', 'e'},
            'b': {'a', 'c', 'd'},
            'c': {'b', 'd'},
            'd': {'b', 'c'},
            'e': {'a'},
        }, {
            'A': {'B'},
            'B': {'A', 'C', 'D', 'E'},
            'C': {'B', 'D', 'E'},
            'D': {'B', 'C'},
            'E': {'B', 'C'},
        }, s3_threshold=None) # set threshold to None so all moves are valid
        
    def test_s3_initial(self):
        self.assertEqual(self.s3_sagrow._s3_frac, [0, 0])

    def test_add_s3_update(self):
        self.s3_sagrow._make_move(0)
        self.assertEqual(self.s3_sagrow._s3_frac, [1, 1])
        self.s3_sagrow._make_move(1)
        self.assertEqual(self.s3_sagrow._s3_frac, [2, 2])
        self.s3_sagrow._make_move(2)
        self.assertEqual(self.s3_sagrow._s3_frac, [4, 4])
        self.s3_sagrow._make_move(3)
        self.assertEqual(self.s3_sagrow._s3_frac, [4, 7])

    def test_add_multiple_s3_update(self):
        self.s3_sagrow._make_move(1)
        self.s3_sagrow._make_move(2)
        self.s3_sagrow._make_move(3)
        self.assertEqual(self.s3_sagrow._s3_frac, [1, 2])
        self.s3_sagrow._make_move(0)
        self.assertEqual(self.s3_sagrow._s3_frac, [4, 7])

    def test_remove_s3_update(self):
        self.s3_sagrow._make_move(0)
        self.s3_sagrow._make_move(1)
        self.s3_sagrow._make_move(2)
        self.s3_sagrow._make_move(3)
        self.assertEqual(self.s3_sagrow._s3_frac, [4, 7])
        self.s3_sagrow._make_move(2)
        self.assertEqual(self.s3_sagrow._s3_frac, [2, 5])
        self.s3_sagrow._make_move(3)
        self.assertEqual(self.s3_sagrow._s3_frac, [2, 2])
        self.s3_sagrow._make_move(1)
        self.assertEqual(self.s3_sagrow._s3_frac, [1, 1])

    def test_can_remove_last_block(self):
        self.s3_sagrow._make_move(0)
        self.assertEqual(self.s3_sagrow._s3_frac, [1, 1])
        self.assertTrue(self.s3_sagrow._move_is_valid(0))

    def test_remove_multiple_s3_update(self):
        self.s3_sagrow._make_move(0)
        self.s3_sagrow._make_move(1)
        self.s3_sagrow._make_move(2)
        self.s3_sagrow._make_move(3)
        self.assertEqual(self.s3_sagrow._s3_frac, [4, 7])
        self.s3_sagrow._make_move(0)
        self.assertEqual(self.s3_sagrow._s3_frac, [1, 2])

    def test_double_add_s3_update(self):
        self.s3_sagrow._make_move(0)
        self.assertEqual(self.s3_sagrow._s3_frac, [1, 1])
        self.s3_sagrow._make_move(3)
        self.assertEqual(self.s3_sagrow._s3_frac, [1, 3])
        # adding a node pair that already exists
        self.s3_sagrow._make_move(4)
        self.assertEqual(self.s3_sagrow._s3_frac, [2, 5])

    def test_double_remove_s3_update(self):
        # when removing a node pair that's covered by another block
        self.s3_sagrow._make_move(1)
        self.s3_sagrow._make_move(2)
        self.s3_sagrow._make_move(3)
        self.s3_sagrow._make_move(4)
        self.assertEqual(self.s3_sagrow._s3_frac, [1, 3])
        self.s3_sagrow._make_move(4)
        self.assertEqual(self.s3_sagrow._s3_frac, [1, 2])

    def test_s3_frac_denom_zero(self):
        self.s3_sagrow._s3_threshold = None
        self.assertTrue(self.s3_sagrow._move_is_valid(2))
        self.s3_sagrow._make_move(2)
        self.assertTrue(self.s3_sagrow._move_is_valid(3))
        self.s3_sagrow._s3_threshold = 0
        self.assertFalse(self.s3_sagrow._move_is_valid(3))
        
    def test_s3_threshold_adding(self):
        self.s3_sagrow._s3_threshold = 1.0
        self.assertTrue(self.s3_sagrow._move_is_valid(0))
        self.s3_sagrow._make_move(0)
        self.assertTrue(self.s3_sagrow._move_is_valid(1))
        self.s3_sagrow._make_move(1)
        self.assertTrue(self.s3_sagrow._move_is_valid(2))
        self.s3_sagrow._make_move(2)
        self.assertFalse(self.s3_sagrow._move_is_valid(3))
        self.s3_sagrow._s3_threshold = 0.55
        self.assertTrue(self.s3_sagrow._move_is_valid(3))
        self.s3_sagrow._s3_threshold = 0.6
        self.assertFalse(self.s3_sagrow._move_is_valid(3))

    def get_blank_sagrow(self):
        return SimAnnealGrow(list(), dict(), dict(), s3_threshold=None)
        
    def test_standard_p_func(self):
        self.assertEqual(self.blank_sagrow._p_func, STANDARD_P_FUNC)
        self.assertEqual(self.blank_sagrow._p(5, 3, 1), 1)
        self.assertEqual(self.blank_sagrow._p(3, 5, 1), math.exp(-2))
        self.assertEqual(self.blank_sagrow._p(5, 3, 0.5), 1)
        self.assertEqual(self.blank_sagrow._p(3, 5, 0.5), math.exp(-4))
        self.assertEqual(self.blank_sagrow._p(5, 3, 0), 1)
        self.assertEqual(self.blank_sagrow._p(3, 5, 0), 0)
        
    def test_always_p_func(self):
        self.blank_sagrow._p_func = ALWAYS_P_FUNC
        self.assertEqual(self.blank_sagrow._p(5, 3, 1), 1)
        self.assertEqual(self.blank_sagrow._p(3, 5, 1), 1)
        self.assertEqual(self.blank_sagrow._p(5, 3, 0.5), 1)
        self.assertEqual(self.blank_sagrow._p(3, 5, 0.5), 1)
        self.assertEqual(self.blank_sagrow._p(5, 3, 0), 1)
        self.assertEqual(self.blank_sagrow._p(3, 5, 0), 1)
        
    def test_energy(self):
        self.assertEqual(energy(8), -8)
    
if __name__ == '__main__':
    unittest.main()
