from tfsm import *
from tfsm_ds import *
import copy
import random
import ast

def ds_rf_uniform_tfsm_derivation(spec_tfsm, impl_tfsm, u, v):
    st_number = len(spec_tfsm.tfsm)
    in_number = 2
    max_len = st_number + max(spec_tfsm.ell, impl_tfsm.ell)
    for l in range(1, max_len):
        i_power_l = in_number ** l
        str_format = '{0:0' + str(l) + 'b}'
        for k in range(0, i_power_l):
            seq = str_format.format(k)
            timed_sequence = []
            for j in range(0, len(seq)):
                i = "i" + str(seq[j])
                theta = float(v - u) / 2
                t = (u + theta) * (j + 1)
                timed_sequence.append((i, t))
            tis = TimedSequence(timed_sequence)
            out_spec = spec_tfsm.generate_output_projections('s0', tis)
            out_impl = impl_tfsm.generate_output_projections('s0', tis)
            if out_spec != out_impl:
                return tis
    return None

class TestSuiteGeneration:
    def __init__(self, tfsm_spec):
        self.tfsm_spec = tfsm_spec
        return

    def derive_output_tran_mutants(self, spec_tran):
        output_mutants = list()
        for o_num in range(0, self.tfsm_spec.outputs_number):
            if o_num != spec_tran.output:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.output = int(o_num)
                output_mutants.append(tfsm_mut)
        return output_mutants

    def derive_transition_tran_mutants(self, spec_tran):
        tran_mutants = list()
        for mut_s in range(0, self.tfsm_spec.states_number):
            if mut_s != spec_tran.end_state:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.end_state = int(mut_s)
                tran_mutants.append(tfsm_mut)
        return tran_mutants

    def derive_delay_tran_mutants(self, spec_tran):
        delay_mutants = list()
        for mut_d in range(1, self.tfsm_spec.d):
            if mut_d != spec_tran.delay:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.delay = int(mut_d)
                delay_mutants.append(tfsm_mut)
        return delay_mutants

    def derive_guard_tran_mutants(self, spec_tran):
        guard_mutants = list()
        for o_num in range(0, self.tfsm_spec.outputs_number):
            if o_num != spec_tran.output:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.output = int(o_num)
                guard_mutants.append(tfsm_mut)
        for mut_s in range(0, self.tfsm_spec.states_number):
            if mut_s != spec_tran.end_state:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.end_state = int(mut_s)
                guard_mutants.append(tfsm_mut)
        for mut_d in range(1, self.tfsm_spec.d):
            if mut_d != spec_tran.delay:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.delay = int(mut_d)
                guard_mutants.append(tfsm_mut)
        return guard_mutants

    def derive_fo_mutants(self, mut_type):
        fo_mutants = list()
        for tran_name in self.tfsm_spec.tfsm.transition_dict.keys():
            spec_tran = self.tfsm_spec.tfsm.transition_dict[tran_name]
            tran_mutants = list()
            if mut_type == "output":
                tran_mutants = self.derive_output_tran_mutants(spec_tran)
            elif mut_type == "transition":
                tran_mutants = self.derive_transition_tran_mutants(spec_tran)
            elif mut_type == "delay":
                tran_mutants = self.derive_delay_tran_mutants(spec_tran)
            elif mut_type == "guard":
                for p in range(self.tfsm_spec.u+1, self.tfsm_spec.v):
                    split_tran1 = copy.deepcopy(spec_tran)
                    split_tran1.transition_name = 'e_1_' + str(spec_tran.transition_name[1:])
                    split_tran1.time_guard = (self.tfsm_spec.u, p)
                    split_tran2 = copy.deepcopy(spec_tran)
                    split_tran2.transition_name = 'e_2_' + str(spec_tran.transition_name[1:])
                    split_tran2.time_guard = (p, self.tfsm_spec.v)
                    spec_split_tfsm = copy.deepcopy(self.tfsm_spec)
                    spec_split_tfsm.tfsm.transition_dict.pop(spec_tran.transition_name)
                    spec_split_tfsm.tfsm.tfsm[spec_tran.start_state][spec_tran.input].pop(spec_tran.time_guard)
                    spec_split_tfsm.tfsm.transition_dict[split_tran1.transition_name] = split_tran1
                    spec_split_tfsm.tfsm.tfsm[split_tran1.start_state][split_tran1.input][split_tran1.time_guard] = split_tran1
                    spec_split_tfsm.tfsm.transition_dict[split_tran2.transition_name] = split_tran2
                    spec_split_tfsm.tfsm.tfsm[split_tran2.start_state][split_tran2.input][split_tran2.time_guard] = split_tran2
                    tsg_split_tfsm = TestSuiteGeneration(spec_split_tfsm)
                    tran_guard_mutants = tsg_split_tfsm.derive_guard_tran_mutants(split_tran1)
                    fo_mutants += tran_guard_mutants
                    tran_guard_mutants = tsg_split_tfsm.derive_guard_tran_mutants(split_tran2)
                    fo_mutants += tran_guard_mutants
            fo_mutants += tran_mutants
        return fo_mutants



class TestSuiteGeneration_RaceFree:
    def __init__(self, tfsm_spec):
        self.tfsm_spec = tfsm_spec
        return

    def extend_untimed_test_suite(self, testSuite):
        testSuite_mean = list()
        testSuite_random = list()
        for test_seq in testSuite:
            ttt_mean_seq = TimedSequence([])
            ttt_random_seq = TimedSequence([])
            n_s = self.tfsm_spec.initial_state
            for i_num in test_seq:
                input = "i" + str(i_num)
                for timed_guard in self.tfsm_spec.tfsm[n_s][input]:
                    tran = self.tfsm_spec.tfsm[n_s][input][timed_guard]
                    timestamp_mean = (float(tran.time_guard[0]) + float(tran.time_guard[1])) / 2
                    timestamp_random = random.uniform(tran.time_guard[0], tran.time_guard[1])
                    ttt_mean_seq.append_to_timed_sequence(TimedSequence([(tran.input, timestamp_mean)]))
                    ttt_random_seq.append_to_timed_sequence(TimedSequence([(tran.input, timestamp_random)]))
            testSuite_mean.append(ttt_mean_seq)
            testSuite_random.append(ttt_random_seq)
        return testSuite_mean, testSuite_random

    def is_killed_by_ts(self, impl_tfsm, ts):
        for tis in ts:
            out_spec = self.tfsm_spec.generate_output_projections('s0', tis)
            out_impl = impl_tfsm.generate_output_projections('s0', tis)
            if out_spec != out_impl:
                return True
        return False

    def mean_extension_for_untimed_test_suite(self, test_suite, is_random=False):
        timed_test_suite = []
        for test_seq in test_suite:
            timed_test_suite.append(self.mean_extension_for_untimed_sequence(test_seq, is_random))
        return timed_test_suite

    def mean_extension_for_untimed_sequence(self, input_seq, is_random=False):
        timed_input_seq = TimedSequence([])
        n_s = self.tfsm_spec.initial_state
        for i in input_seq:
            for timed_guard in self.tfsm_spec.tfsm.tfsm[n_s][i]:
                tran = self.tfsm_spec.tfsm.tfsm[n_s][i][timed_guard]
                if is_random:
                    timestamp_mean = random.uniform(tran.time_guard[0], tran.time_guard[1])
                else:
                    timestamp_mean = (float(tran.time_guard[0]) + float(tran.time_guard[1])) / 2
                timed_input_seq.append_to_timed_sequence(TimedSequence([(i, timestamp_mean)]))
                n_s = int(tran.end_state)
        return timed_input_seq

    def derive_output_tran_mutants(self, spec_tran):
        output_mutants = list()
        output_test_sequences = list()
        for o_num in range(0, self.tfsm_spec.outputs_number):
            if o_num != spec_tran.output:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.output = int(o_num)
                output_mutants.append(tfsm_mut)
                ts = list(self.tfsm_spec.fsm.transfer_sequences[mut_tran.start_state]) + [mut_tran.input]
                output_test_sequences.append(self.mean_extension_for_untimed_sequence(ts))
        return output_mutants, output_test_sequences

    def derive_transition_tran_mutants(self, spec_tran):
        tran_mutants = list()
        tran_test_sequences = list()
        for mut_s in range(0, self.tfsm_spec.states_number):
            if mut_s != spec_tran.end_state:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.end_state = int(mut_s)
                #if tfsm_mut.tfsm.is_race_free():
                tfsm_mut.generate_reverse_transition_dict()
                if tfsm_mut.is_race_free_tran(mut_tran):
                    max_len = self.tfsm_spec.states_number + self.tfsm_spec.tfsm.ell
                    rf_tree = RaceFree_TruncatedTree_for_testing(self.tfsm_spec, tfsm_mut, 0, 0, max_len)
                    tis_ds = rf_tree.derive_bfs_tree()
                    if tis_ds:
                        tran_mutants.append(tfsm_mut)
                        tran_test_sequences.append(tis_ds)
        return tran_mutants, tran_test_sequences

    def derive_delay_tran_mutants(self, spec_tran):
        delay_mutants = list()
        delay_test_sequences = list()
        for mut_d in range(1, self.tfsm_spec.d):
            if mut_d != spec_tran.delay:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.delay = int(mut_d)
                #if tfsm_mut.tfsm.is_race_free():
                tfsm_mut.generate_reverse_transition_dict()
                if tfsm_mut.is_race_free_tran(mut_tran):
                    max_len = self.tfsm_spec.states_number
                    rf_tree = RaceFree_TruncatedTree_for_testing(self.tfsm_spec, tfsm_mut, spec_tran.start_state, mut_tran.start_state, max_len)
                    tis_ds = rf_tree.derive_bfs_tree()
                    if tis_ds:
                        ts_proj = list(self.tfsm_spec.fsm.transfer_sequences[mut_tran.start_state])
                        ts = self.mean_extension_for_untimed_sequence(ts_proj)
                        ts.append_to_timed_sequence(tis_ds)
                        spec_out = self.tfsm_spec.derive_output_sequence(0, ts)
                        mut_out = tfsm_mut.derive_output_sequence(0, ts)
                        if spec_out != mut_out:
                            delay_mutants.append(tfsm_mut)
                            delay_test_sequences.append(ts)
        return delay_mutants, delay_test_sequences

    def derive_guard_tran_mutants(self, spec_tran):
        guard_mutants = list()
        guard_test_sequences = list()
        max_len = self.tfsm_spec.states_number + self.tfsm_spec.tfsm.ell - 1
        for o_num in range(0, self.tfsm_spec.outputs_number):
            if o_num != spec_tran.output:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.output = int(o_num)
                rf_tree = RaceFree_TruncatedTree_for_testing(self.tfsm_spec, tfsm_mut, 0, 0, max_len)
                tis_ds = rf_tree.derive_bfs_tree()
                if tis_ds:
                    guard_mutants.append(tfsm_mut)
                    guard_test_sequences.append(tis_ds)
        for mut_s in range(0, self.tfsm_spec.states_number):
            if mut_s != spec_tran.end_state:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.end_state = int(mut_s)
                tfsm_mut.generate_reverse_transition_dict()
                if tfsm_mut.is_race_free_tran(mut_tran):
                    rf_tree = RaceFree_TruncatedTree_for_testing(self.tfsm_spec, tfsm_mut, 0, 0, max_len)
                    tis_ds = rf_tree.derive_bfs_tree()
                    if tis_ds:
                        guard_mutants.append(tfsm_mut)
                        guard_test_sequences.append(tis_ds)
        for mut_d in range(1, self.tfsm_spec.d):
            if mut_d != spec_tran.delay:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.delay = int(mut_d)
                tfsm_mut.generate_reverse_transition_dict()
                if tfsm_mut.is_race_free_tran(mut_tran):
                    max_len = self.tfsm_spec.states_number
                    rf_tree = RaceFree_TruncatedTree_for_testing(self.tfsm_spec, tfsm_mut, 0, 0, max_len)
                    tis_ds = rf_tree.derive_bfs_tree()
                    if tis_ds:
                        guard_mutants.append(tfsm_mut)
                        guard_test_sequences.append(tis_ds)
        return guard_mutants, guard_test_sequences

    def derive_fo_mutants(self, mut_type):
        fo_mutants = list()
        fo_tests = list()
        for tran_name in self.tfsm_spec.tfsm.transition_dict.keys():
            print("tran_name =", tran_name)
            spec_tran = self.tfsm_spec.tfsm.transition_dict[tran_name]
            tran_mutants = list()
            tran_tests = list()
            if mut_type == "output":
                tran_mutants, tran_tests = self.derive_output_tran_mutants(spec_tran)
            elif mut_type == "transition":
                tran_mutants, tran_tests = self.derive_transition_tran_mutants(spec_tran)
            elif mut_type == "delay":
                tran_mutants, tran_tests = self.derive_delay_tran_mutants(spec_tran)
            elif mut_type == "guard":
                #for p in range(self.tfsm_spec.u+1, self.tfsm_spec.v):
                for p in range(spec_tran.time_guard[0] + 1, spec_tran.time_guard[1]):
                    split_tran1 = copy.deepcopy(spec_tran)
                    split_tran1.transition_name = 'e_1_' + str(spec_tran.transition_name[1:])
                    split_tran1.time_guard = (spec_tran.time_guard[0], p)
                    split_tran2 = copy.deepcopy(spec_tran)
                    split_tran2.transition_name = 'e_2_' + str(spec_tran.transition_name[1:])
                    split_tran2.time_guard = (p, spec_tran.time_guard[1])
                    spec_split_tfsm = copy.deepcopy(self.tfsm_spec)
                    spec_split_tfsm.tfsm.transition_dict.pop(spec_tran.transition_name)
                    spec_split_tfsm.tfsm.tfsm[spec_tran.start_state][spec_tran.input].pop(spec_tran.time_guard)
                    spec_split_tfsm.tfsm.transition_dict[split_tran1.transition_name] = split_tran1
                    spec_split_tfsm.tfsm.tfsm[split_tran1.start_state][split_tran1.input][split_tran1.time_guard] = split_tran1
                    spec_split_tfsm.tfsm.transition_dict[split_tran2.transition_name] = split_tran2
                    spec_split_tfsm.tfsm.tfsm[split_tran2.start_state][split_tran2.input][split_tran2.time_guard] = split_tran2
                    tsg_split_tfsm = TestSuiteGeneration_RaceFree(spec_split_tfsm)
                    tran_guard_mutants, tran_guard_tests = tsg_split_tfsm.derive_guard_tran_mutants(split_tran1)
                    fo_mutants += tran_guard_mutants
                    tran_tests += tran_guard_tests
                    tran_guard_mutants, tran_guard_tests = tsg_split_tfsm.derive_guard_tran_mutants(split_tran2)
                    fo_mutants += tran_guard_mutants
                    tran_tests += tran_guard_tests
            fo_mutants += tran_mutants
            fo_tests += tran_tests
        return fo_mutants, fo_tests

    def derive_fo_output_mutants(self):
        output_mutants = list()
        output_test_sequences = list()
        for tran in self.tfsm_spec.tfsm.transition_dict:
            tfsm_mut = copy.deepcopy(self.tfsm_spec)
            mut_tran = tfsm_mut.tfsm.transition_dict[tran]
            mut_tran.output = (mut_tran.output + 1) % 2
            output_mutants.append(tfsm_mut)
            ts = list(self.tfsm_spec.fsm.transfer_sequences[mut_tran.start_state]) + [mut_tran.input]
            output_test_sequences.append(self.mean_extension_for_untimed_sequence(ts))
        return output_mutants, output_test_sequences

    def derive_fo_tran_mutants(self):
        tran_mutants = list()
        tran_test_sequences = list()
        for tran in self.tfsm_spec.tfsm.transition_dict:
            spec_tran = self.tfsm_spec.tfsm.transition_dict[tran]
            for mut_s in range(0, self.tfsm_spec.states_number):
                if mut_s != spec_tran.end_state:
                    tfsm_mut = copy.deepcopy(self.tfsm_spec)
                    mut_tran = tfsm_mut.tfsm.transition_dict[tran]
                    mut_tran.end_state = int(mut_s)
                    if tfsm_mut.tfsm.is_race_free():
                        max_len = self.tfsm_spec.states_number
                        rf_tree = RaceFree_TruncatedTree_for_testing(self.tfsm_spec, tfsm_mut, spec_tran.end_state, mut_tran.end_state, max_len)
                        tis_ds = rf_tree.derive_bfs_tree()
                        if tis_ds:
                            ts_proj = list(self.tfsm_spec.fsm.transfer_sequences[mut_tran.start_state]) + [spec_tran.input]
                            ts = self.mean_extension_for_untimed_sequence(ts_proj)
                            ts.append_to_timed_sequence(tis_ds)
                            tran_mutants.append(tfsm_mut)
                            tran_test_sequences.append(ts)
        return tran_mutants, tran_test_sequences

    def derive_fo_delay_mutants(self):
        delay_mutants = list()
        delay_test_sequences = list()
        for tran in self.tfsm_spec.tfsm.transition_dict:
            spec_tran = self.tfsm_spec.tfsm.transition_dict[tran]
            for mut_d in range(1, self.tfsm_spec.d):
                if mut_d != spec_tran.delay:
                    tfsm_mut = copy.deepcopy(self.tfsm_spec)
                    mut_tran = tfsm_mut.tfsm.transition_dict[tran]
                    mut_tran.delay = int(mut_d)
                    if tfsm_mut.tfsm.is_race_free():
                        max_len = self.tfsm_spec.states_number
                        rf_tree = RaceFree_TruncatedTree_for_testing(self.tfsm_spec, tfsm_mut, spec_tran.start_state, mut_tran.start_state, max_len)
                        tis_ds = rf_tree.derive_bfs_tree()
                        if tis_ds:
                            ts_proj = list(self.tfsm_spec.fsm.transfer_sequences[mut_tran.start_state])
                            ts = self.mean_extension_for_untimed_sequence(ts_proj)
                            ts.append_to_timed_sequence(tis_ds)
                            spec_out = self.tfsm_spec.derive_output_sequence(0, ts)
                            mut_out = tfsm_mut.derive_output_sequence(0, ts)
                            if spec_out != mut_out:
                                delay_mutants.append(tfsm_mut)
                                delay_test_sequences.append(ts)
        return delay_mutants, delay_test_sequences

    def derive_fo_guard_mutants(self):
        grid_tfsm = copy.deepcopy(self.tfsm_spec)
        j = 0
        for tran in self.tfsm_spec.tfsm.transition_dict.keys():
            spec_tran = copy.deepcopy(self.tfsm_spec.tfsm.transition_dict[tran])
            grid_tfsm.tfsm.transition_dict.pop(tran)
            grid_tfsm.tfsm.tfsm[spec_tran.start_state][spec_tran.input].pop(spec_tran.time_guard)
            for p in range(self.tfsm_spec.u+1, self.tfsm_spec.v):
                grid_tran = copy.deepcopy(spec_tran)
                grid_tran.transition_name = "g"+str(j)
                grid_tran.time_guard = (int(p)-1, int(p))
                grid_tfsm.tfsm.transition_dict[grid_tran.transition_name] = grid_tran
                grid_tfsm.tfsm.tfsm[spec_tran.start_state][spec_tran.input][grid_tran.transition_name] = grid_tran
                j += 1
        guard_mutants = list()
        guard_test_sequences = list()
        for tran in self.tfsm_spec.tfsm.transition_dict:
            spec_tran = self.tfsm_spec.tfsm.transition_dict[tran]
            for p in range(spec_tran.time_guard[0]+1,spec_tran.time_guard[1]):
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict.pop(tran)
                tfsm_mut.tfsm.tfsm[spec_tran.start_state][spec_tran.input].pop(spec_tran.time_guard)
                mut_tran1 = copy.deepcopy(mut_tran)
                mut_tran1
                mut_tran2 = copy.deepcopy(mut_tran)
        return

    def generate_test_suites(self, rf_spec_tfsm, type):
        mutants = list()
        testSuite = list()
        for tran in self.tfsm_spec.transition_dict:
            #print("tran =", tran)
            impl_tfsm_list, ds_list = self.derive_first_order_mutant(self.tfsm_spec.transition_dict[tran], rf_spec_tfsm, type)
            for k in range(0, len(impl_tfsm_list)):
                mutants.append(impl_tfsm_list[k])
                if not self.is_killed_by_ts(impl_tfsm_list[k], testSuite):
                    testSuite.append(ds_list[k])
        return mutants, testSuite

    def derive_test_cover(self, mutants_list, test_lists):
        table_test_cover = dict()
        for test_seq in test_lists:
            if not test_seq in table_test_cover:
                table_test_cover[test_seq] = set()
                for index_mut, curr_mut in enumerate(mutants_list):
                    out_spec = self.tfsm_spec.derive_output_sequence(0, test_seq)
                    out_mut = curr_mut.derive_output_sequence(0, test_seq)
                    if out_spec and out_mut and out_spec != out_mut:
                        table_test_cover[test_seq].add(index_mut)
        return table_test_cover

    def derive_complete_test_gradient(self, mutants_list, test_lists):
        table_test_cover = self.derive_test_cover(mutants_list, test_lists)
        uncovered_mutants = set(range(len(mutants_list)))
        complete_test_suite = list()
        while len(uncovered_mutants) > 0:
            test_with_max_cover = max(table_test_cover, key=lambda ts: len(table_test_cover[ts]))
            mutants_to_remove = set(table_test_cover[test_with_max_cover])
            if len(mutants_to_remove) == 0:
                print("STOP")
            complete_test_suite.append(test_with_max_cover)
            uncovered_mutants.difference_update(mutants_to_remove)
            for ts in table_test_cover.keys():
                table_test_cover[ts].difference_update(mutants_to_remove)
        return complete_test_suite

    def apply_ts_to_mutants(self, ts, mutant_list):
        killed_mutants = 0
        for impl_tfsm in mutant_list:
            is_killed = False
            k = 0
            while (k < len(ts)) and (not is_killed):
                tis = ts[k]
                out_spec = self.tfsm_spec.derive_output_sequence(0, tis)
                out_impl = impl_tfsm.derive_output_sequence(0, tis)
                if out_spec != out_impl:
                    is_killed = True
                k += 1
            if is_killed:
                killed_mutants += 1
        return killed_mutants