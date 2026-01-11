import math
from collections import deque
from itertools import permutations
import random
import copy
import ast

class State_pair_fsm:
    def __init__(self, s1, s2):
        self.s1 = int(s1)
        self.s2 = int(s2)
        self.is_distinguished = False
        self.shortest_ds = None
        return

    def __hash__(self):
        return hash((self.s1, self.s2))

class FSM:
    def __init__(self, fsm_file):
        self.fsm_file = str(fsm_file)
        self.transition_dict = {}
        return

    def parse_fsm(self):
        file = open(self.fsm_file, 'r')
        j = 0
        for line in file:
            parts = line.strip().split()
            if not parts or parts[0] in {'F', 'S', 'I', 'O', 'P'}:
                continue
            if parts[0] == 's':
                self.states_number = int(parts[1])
                continue
            if parts[0] == 'i':
                self.inputs_number = int(parts[1])
                continue
            if parts[0] == 'o':
                self.outputs_number = int(parts[1])
                continue
            if parts[0] == 'n0':
                self.initial_state = int(parts[1])
                continue
            if len(parts) >= 4:
                j += 1
                start_state_num, input_num, end_state_num, output_num = parts[:4]
                start_state = int(start_state_num)
                input_d = int(input_num)
                output = int(output_num)
                end_state = int(end_state_num)
                if start_state not in self.transition_dict:
                    self.transition_dict[start_state] = {}
                self.transition_dict[start_state][input_d] = [output, end_state]
        return

    def derive_reverse_fsm(self):
        self.reverse_transition_dict = {(s, i): [] for s in range(self.states_number) for i in range(self.inputs_number)}
        for s in range(0, self.states_number):
            for i in range(0, self.inputs_number):
                [o, n_s] = self.transition_dict[s][i]
                self.reverse_transition_dict[(n_s, i)].append([s, i, o, n_s])
        return

    def bottom_up_alg_step(self, k_dist, k_dist_shortest):
        k_plus_1_dist = set(k_dist)
        k_plus_1_dist_shortest = set()
        for pair in k_dist_shortest:
            s1 = pair[0]
            s2 = pair[1]
            for i in range(0, self.inputs_number):
                for [s1_prev, _, _, _] in self.reverse_transition_dict[(s1, i)]:
                    for [s2_prev, _, _, _] in self.reverse_transition_dict[(s2, i)]:
                        prev_pair = tuple(sorted((s1_prev, s2_prev)))
                        if (s1_prev != s2_prev) and (not prev_pair in k_plus_1_dist):
                            k_plus_1_dist.add(prev_pair)
                            k_plus_1_dist_shortest.add(prev_pair)
                            self.shortest_DSs[prev_pair] = [i] + list(self.shortest_DSs[pair])
        return k_plus_1_dist, k_plus_1_dist_shortest

    def derive_shortest_DSs_bottom_up(self):
        self.shortest_DSs = {(s1, s2) : list() for s1 in range(0, self.states_number) for s2 in range(s1+1, self.states_number)}
        one_dist = set()
        for pair in self.shortest_DSs.keys():
            s1 = pair[0]
            s2 = pair[1]
            is_states_one_dist = False
            i = 0
            while (i < self.inputs_number) and (not is_states_one_dist):
                [o1, _] = self.transition_dict[s1][i]
                [o2, _] = self.transition_dict[s2][i]
                if o1 != o2:
                    is_states_one_dist = True
                    self.shortest_DSs[(s1, s2)] = [i]
                    one_dist.add((s1, s2))
                i += 1
        if len(one_dist) < (self.states_number * (self.states_number - 1)) / 2:
            k_dist = set(one_dist)
            k_dist_shortest = set(one_dist)
            stop_flag = False
            while not stop_flag:
                k_plus_1_dist, k_plus_1_dist_shortest = self.bottom_up_alg_step(k_dist, k_dist_shortest)
                if len(k_plus_1_dist) == len(k_dist):
                    stop_flag = True
                else:
                    k_dist = set(k_plus_1_dist)
                    k_dist_shortest = set(k_plus_1_dist_shortest)
        return

    def find_reachable_seq(self):
        self.transfer_sequences = dict()
        for s in self.transition_dict:
            self.transfer_sequences[s] = list()
        visited_states = [0] * len(self.transfer_sequences)
        init_state = int(self.initial_state)
        queue = [init_state]
        visited_states[init_state] = 1
        while len(queue) > 0:
            curr_state = queue.pop(0)
            for i in self.transition_dict[curr_state]:
                next_state = int(self.transition_dict[curr_state][i][1])
                if not visited_states[next_state]:
                    visited_states[next_state] = 1
                    self.transfer_sequences[next_state] = list(self.transfer_sequences[curr_state]) + [int(i)]
                    queue.append(next_state)
        return

    def find_ds_for_two_fsms(self, s, fsm2, q):
        N = max(self.states_number, fsm2.states_number)
        in_number = 2
        for l in range(1, N+1):
            i_power_l = in_number ** l
            str_format = '{0:0' + str(l) + 'b}'
            for k in range(0, i_power_l):
                seq = str_format.format(k)
                n_s = int(s)
                n_q = int(q)
                for i in range(0, l-1):
                    n_s = int(self.transition_dict[n_s][int(seq[i])][1])
                    n_q = int(fsm2.transition_dict[n_q][int(seq[i])][1])
                if self.transition_dict[n_s][int(seq[-1])][0] != self.transition_dict[n_q][int(seq[-1])][0]:
                    return list(map(int, seq))
        return None

    def derive_output_sequence(self, input_seq):
        out = list()
        n_s = int(self.initial_state)
        for i in input_seq:
            [o, n_s] = self.transition_dict[n_s][i]
            out.append(o)
        return out

    def find_ds(self, s1, s2):
        st_number = len(self.transition_dict)
        # in_number = len(self.transition_dict["s0"])
        in_number = 2
        for l in range(1, st_number):
            i_power_l = in_number ** l
            str_format = '{0:0' + str(l) + 'b}'
            for k in range(0, i_power_l):
                seq = str_format.format(k)
                n_s1 = int(s1)
                n_s2 = int(s2)
                for i in range(0, l-1):
                    n_s1 = int(self.transition_dict[n_s1][int(seq[i])][1])
                    n_s2 = int(self.transition_dict[n_s2][int(seq[i])][1])
                if self.transition_dict[n_s1][int(seq[-1])][0] != self.transition_dict[n_s2][int(seq[-1])][0]:
                    return seq
        return None

class FSM_FirstOrderTestSuiteGeneration:
    def __init__(self, fsm_spec):
        self.fsm_spec = fsm_spec
        return

    def is_mutant_killed_by_test_suite(self, impl_fsm, testSuite):
        for test_seq in testSuite:
            out_spec = self.fsm_spec.derive_output_sequence(test_seq)
            out_impl = impl_fsm.derive_output_sequence(test_seq)
            if out_spec != out_impl:
                return True
        return False

    def is_mutant_killed_by_test_seq(self, impl_fsm, test_seq):
        out_spec = self.fsm_spec.derive_output_sequence(test_seq)
        out_impl = impl_fsm.derive_output_sequence(test_seq)
        if out_spec != out_impl:
            return True
        return False

    def derive_first_order_output_mutants(self):
        output_mutants = list()
        output_test_sequences = list()
        for s in range(0, self.fsm_spec.states_number):
            for i in range(0, self.fsm_spec.inputs_number):
                [o, n_s] = self.fsm_spec.transition_dict[s][i]
                o_mut = (o + 1) % 2
                fsm_impl = copy.deepcopy(self.fsm_spec)
                fsm_impl.transition_dict[s][i] = [o_mut, n_s]
                output_mutants.append(fsm_impl)
                ts = list(self.fsm_spec.transfer_sequences[s]) + [i]
                output_test_sequences.append(ts)
        return output_mutants, output_test_sequences

    def derive_first_order_transition_mutants(self):
        transition_mutants = list()
        transition_test_sequences = list()
        for s in range(0, self.fsm_spec.states_number):
            for i in range(0, self.fsm_spec.inputs_number):
                [o, n_s] = self.fsm_spec.transition_dict[s][i]
                for mut_n_s in range(0, self.fsm_spec.states_number):
                    if mut_n_s != n_s:
                        fsm_impl = copy.deepcopy(self.fsm_spec)
                        fsm_impl.transition_dict[s][i] = [o, mut_n_s]
                        pair = tuple(sorted((n_s, mut_n_s)))
                        if self.fsm_spec.shortest_DSs[pair]:
                            transition_mutants.append(fsm_impl)
                            ts = list(self.fsm_spec.transfer_sequences[s]) + [i] + list(self.fsm_spec.shortest_DSs[pair])
                            transition_test_sequences.append(ts)
        return transition_mutants, transition_test_sequences

    def derive_test_cover(self, mutants_list, test_lists):
        table_test_cover = dict()
        for test_seq in test_lists:
            if not str(test_seq) in table_test_cover:
                table_test_cover[str(test_seq)] = set()
                for index_mut, curr_mut in enumerate(mutants_list):
                    if self.is_mutant_killed_by_test_seq(curr_mut, test_seq):
                        table_test_cover[str(test_seq)].add(index_mut)
        return table_test_cover

    def derive_complete_test_gradient(self, mutants_list, test_lists):
        table_test_cover = self.derive_test_cover(mutants_list, test_lists)
        uncovered_mutants = set(range(len(mutants_list)))
        complete_test_suite = list()
        while len(uncovered_mutants) > 0:
            test_with_max_cover = max(table_test_cover, key=lambda ts: len(table_test_cover[ts]))
            mutants_to_remove = set(table_test_cover[test_with_max_cover])
            complete_test_suite.append(ast.literal_eval(test_with_max_cover))
            uncovered_mutants.difference_update(mutants_to_remove)
            for ts in table_test_cover.keys():
                table_test_cover[ts].difference_update(mutants_to_remove)
        return complete_test_suite

    def apply_ts_to_mutants(self, ts, mutant_list):
        killed_mutants = 0
        for impl_fsm in mutant_list:
            is_killed = False
            k = 0
            while (k < len(ts)) and (not is_killed):
                tis = ts[k]
                out_spec = self.fsm_spec.derive_output_sequence(tis)
                out_impl = impl_fsm.derive_output_sequence(tis)
                if out_spec != out_impl:
                    is_killed = True
                k += 1
            if is_killed:
                killed_mutants += 1
        return float(killed_mutants) / float(len(mutant_list))