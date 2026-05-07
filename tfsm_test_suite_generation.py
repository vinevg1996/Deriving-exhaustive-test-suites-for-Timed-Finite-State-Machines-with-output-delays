from tfsm import *
from tfsm_ds import *
import copy
import random
import ast

class TestSuiteGeneration:
    def __init__(self, tfsm_spec):
        self.tfsm_spec = tfsm_spec
        self.fo_mutants_number = 0
        return

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
        if len(mutant_list):
            return float(killed_mutants) / float(len(mutant_list))
        else:
            return 1.0

    def transition_extend_to_timed_test_suite(self, test_suite, strategy):
        timed_test_suite = []
        for test_seq in test_suite:
            epsilon = 1 / float(len(test_seq))
            timed_test_seq = TimedSequence([])
            for tran_name in test_seq:
                tran = self.tfsm_spec.tfsm.transition_dict[tran_name]
                if strategy == "left":
                    timestamp = float(tran.time_guard[0])
                elif strategy == "right":
                    timestamp = float(tran.time_guard[1]) - epsilon
                elif strategy == "mean":
                    timestamp = float(tran.time_guard[0] + tran.time_guard[1]) / 2
                else:
                    timestamp = random.uniform(tran.time_guard[0], tran.time_guard[1])
                timed_test_seq.append_to_timed_sequence(TimedSequence([(tran.input, timestamp)]))
            timed_test_suite.append(timed_test_seq)
        return timed_test_suite

    def extend_to_timed_test_suite(self, test_suite, strategy):
        timed_test_suite = []
        for test_seq in test_suite:
            timed_test_suite.append(self.extend_to_timed_test_sequence(test_seq, strategy))
        return timed_test_suite

    def extend_to_timed_test_sequence(self, input_seq, strategy):
        timed_input_seq = TimedSequence([])
        n_s = self.tfsm_spec.initial_state
        if len(input_seq)==0:
            return TimedSequence([])
        epsilion = 1 / float(2*len(input_seq))
        for i in input_seq:
            for timed_guard in self.tfsm_spec.tfsm.tfsm[n_s][i]:
                tran = self.tfsm_spec.tfsm.tfsm[n_s][i][timed_guard]
                if strategy == "random":
                    timestamp = random.uniform(tran.time_guard[0], tran.time_guard[1])
                elif strategy == "mean":
                    timestamp = (float(tran.time_guard[0]) + float(tran.time_guard[1])) / 2
                elif strategy == "right":
                    timestamp = float(tran.time_guard[1]) - epsilion
                elif strategy == "left":
                    timestamp = float(tran.time_guard[0])
                else:
                    timestamp = -1.0
                timed_input_seq.append_to_timed_sequence(TimedSequence([(i, timestamp)]))
                n_s = int(tran.end_state)
                break
        return timed_input_seq

    def derive_rand_timed_mutant_fast(self, order):
        curr_tfsm_mut = copy.deepcopy(self.tfsm_spec)
        potential_tcs = list()
        for curr_order in range(0, order):
            tran_name = random.choice(list(curr_tfsm_mut.tfsm.transition_dict.keys()))
            spec_tran = curr_tfsm_mut.tfsm.transition_dict[tran_name]
            if (spec_tran.time_guard[1] - spec_tran.time_guard[0]) > 1:
                p = random.randint(spec_tran.time_guard[0] + 1, spec_tran.time_guard[1] - 1)
                split_tran1 = copy.deepcopy(spec_tran)
                split_tran1.transition_name = str(spec_tran.transition_name) + '_1'
                split_tran1.time_guard = (spec_tran.time_guard[0], p)
                split_tran2 = copy.deepcopy(spec_tran)
                split_tran2.transition_name = str(spec_tran.transition_name) + '_2'
                split_tran2.time_guard = (p, spec_tran.time_guard[1])
                curr_tfsm_mut.tfsm.transition_dict.pop(spec_tran.transition_name)
                curr_tfsm_mut.tfsm.tfsm[spec_tran.start_state][spec_tran.input].pop(spec_tran.time_guard)
                curr_tfsm_mut.tfsm.transition_dict[split_tran1.transition_name] = split_tran1
                curr_tfsm_mut.tfsm.tfsm[split_tran1.start_state][split_tran1.input][
                    split_tran1.time_guard] = split_tran1
                curr_tfsm_mut.tfsm.transition_dict[split_tran2.transition_name] = split_tran2
                curr_tfsm_mut.tfsm.tfsm[split_tran2.start_state][split_tran2.input][
                    split_tran2.time_guard] = split_tran2
                choice_tran = random.choice([1, 2])
                if choice_tran == 1:
                    spec_tran = split_tran1
                else:
                    spec_tran = split_tran2
            choice_mut = random.choice(["output", "tran", "delay"])
            if choice_mut == "output":
                curr_o = int(spec_tran.output)
                mut_o = random.randint(0, self.tfsm_spec.outputs_number - 1)
                while mut_o == spec_tran.output:
                    mut_o = random.randint(0, self.tfsm_spec.outputs_number - 1)
                spec_tran.output = int(mut_o)
                untimed_seq = list(self.tfsm_spec.fsm.transfer_sequences[spec_tran.start_state])
                tis_ds = self.extend_to_timed_test_sequence(untimed_seq, "left")
                tis_ds.append_to_timed_sequence(TimedSequence([(spec_tran.input, spec_tran.time_guard[0])]))
                potential_tcs.append(tis_ds)
            elif choice_mut == "tran":
                curr_end_state = int(spec_tran.end_state)
                mut_s = random.randint(0, self.tfsm_spec.states_number - 1)
                while mut_s == spec_tran.end_state:
                    mut_s = random.randint(0, self.tfsm_spec.states_number - 1)
                spec_tran.end_state = int(mut_s)
                tis_ds = self.extend_to_timed_test_sequence(self.tfsm_spec.fsm.transfer_sequences[spec_tran.start_state], "left")
                tis_ds.append_to_timed_sequence(TimedSequence([(spec_tran.input, spec_tran.time_guard[0])]))
                if curr_end_state < mut_s:
                    untimed_tail = list(self.tfsm_spec.fsm.shortest_DSs[curr_end_state, mut_s])
                else:
                    untimed_tail = list(self.tfsm_spec.fsm.shortest_DSs[mut_s, curr_end_state])
                timed_tail = self.extend_to_timed_test_sequence(untimed_tail, "left")
                tis_ds.append_to_timed_sequence(timed_tail)
                potential_tcs.append(tis_ds)
            elif choice_mut == "delay":
                curr_d = int(spec_tran.delay)
                mut_d = random.randint(1, self.tfsm_spec.d)
                while mut_d == spec_tran.delay:
                    mut_d = random.randint(1, self.tfsm_spec.d)
                spec_tran.delay = int(mut_d)
                rf_tree = RaceFree_TruncatedTree_for_testing(self.tfsm_spec, curr_tfsm_mut, spec_tran.start_state, spec_tran.start_state, self.tfsm_spec.tfsm.ell)
                tail_ds = rf_tree.derive_bfs_tree()
                if tail_ds:
                    tis_ds = self.extend_to_timed_test_sequence(self.tfsm_spec.fsm.transfer_sequences[spec_tran.start_state], "left")
                    tis_ds.append_to_timed_sequence(tail_ds)
                    potential_tcs.append(tis_ds)
        for curr_ds in potential_tcs:
            out_spec = self.tfsm_spec.derive_output_sequence(0, curr_ds)
            out_mut = curr_tfsm_mut.derive_output_sequence(0, curr_ds)
            if out_spec and out_mut and out_spec != out_mut:
                return curr_tfsm_mut
        return None

    def derive_rand_timed_mutant(self, order, mutants_type):
        curr_tfsm_mut = copy.deepcopy(self.tfsm_spec)
        for it in range(0, order):
            if mutants_type == "delay":
                self.derive_rand_delay_mutant(curr_tfsm_mut)
            else:
                self.derive_rand_guard_mutant(curr_tfsm_mut)
        max_len = self.tfsm_spec.states_number + self.tfsm_spec.tfsm.ell
        rf_tree = RaceFree_TruncatedTree_for_testing(self.tfsm_spec, curr_tfsm_mut, 0, 0, max_len)
        tis_ds = rf_tree.derive_bfs_tree()
        if tis_ds:
            return curr_tfsm_mut, tis_ds
        else:
            return None, None

    def derive_rand_delay_mutant(self, curr_tfsm_mut):
        tran_name = random.choice(list(curr_tfsm_mut.tfsm.transition_dict.keys()))
        spec_tran = curr_tfsm_mut.tfsm.transition_dict[tran_name]
        mut_d = random.randint(1, curr_tfsm_mut.d)
        while mut_d == spec_tran.delay:
            mut_d = random.randint(1, curr_tfsm_mut.d)
        spec_tran.delay = int(mut_d)
        return

    def derive_rand_guard_mutant(self, curr_tfsm_mut):
        tran_name = random.choice(list(curr_tfsm_mut.tfsm.transition_dict.keys()))
        spec_tran = curr_tfsm_mut.tfsm.transition_dict[tran_name]
        if (spec_tran.time_guard[1] - spec_tran.time_guard[0]) > 1:
            p = random.randint(spec_tran.time_guard[0]+1, spec_tran.time_guard[1]-1)
            split_tran1 = copy.deepcopy(spec_tran)
            split_tran1.transition_name = str(spec_tran.transition_name) + '_1'
            split_tran1.time_guard = (spec_tran.time_guard[0], p)
            split_tran2 = copy.deepcopy(spec_tran)
            split_tran2.transition_name = str(spec_tran.transition_name) + '_2'
            split_tran2.time_guard = (p, spec_tran.time_guard[1])
            curr_tfsm_mut.tfsm.transition_dict.pop(spec_tran.transition_name)
            curr_tfsm_mut.tfsm.tfsm[spec_tran.start_state][spec_tran.input].pop(spec_tran.time_guard)
            curr_tfsm_mut.tfsm.transition_dict[split_tran1.transition_name] = split_tran1
            curr_tfsm_mut.tfsm.tfsm[split_tran1.start_state][split_tran1.input][split_tran1.time_guard] = split_tran1
            curr_tfsm_mut.tfsm.transition_dict[split_tran2.transition_name] = split_tran2
            curr_tfsm_mut.tfsm.tfsm[split_tran2.start_state][split_tran2.input][split_tran2.time_guard] = split_tran2
            choice_tran = random.choice([1, 2])
            if choice_tran == 1:
                spec_tran = split_tran1
            else:
                spec_tran = split_tran2
        choice_mut = random.choice(["output", "tran", "delay"])
        if choice_mut == "output":
            mut_o = random.randint(0, self.tfsm_spec.outputs_number-1)
            while mut_o == spec_tran.output:
                mut_o = random.randint(0, self.tfsm_spec.outputs_number-1)
            spec_tran.output = int(mut_o)
        elif choice_mut == "tran":
            mut_s = random.randint(0, self.tfsm_spec.states_number-1)
            while mut_s == spec_tran.end_state:
                mut_s = random.randint(0, self.tfsm_spec.states_number-1)
            spec_tran.end_state = int(mut_s)
        elif choice_mut == "delay":
            mut_d = random.randint(1, self.tfsm_spec.d)
            while mut_d == spec_tran.delay:
                mut_d = random.randint(1, self.tfsm_spec.d)
            spec_tran.delay = int(mut_d)
        return

    def derive_all_fo_mutants(self, mut_type):
        fo_mutants = list()
        for tran_name in self.tfsm_spec.tfsm.transition_dict.keys():
            #print("tran_name =", tran_name)
            spec_tran = self.tfsm_spec.tfsm.transition_dict[tran_name]
            tran_mutants = list()
            if mut_type == "output":
                tran_mutants = self.derive_all_output_tran_mutants(spec_tran)
            elif mut_type == "transition":
                tran_mutants = self.derive_all_transition_tran_mutants(spec_tran)
            elif mut_type == "delay":
                tran_mutants = self.derive_all_delay_tran_mutants(spec_tran)
            elif mut_type == "guard":
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
                    if spec_split_tfsm.tfsm.is_race_free():
                        tsg_split_tfsm = TestSuiteGeneration(spec_split_tfsm)
                        tran_guard_mutants = tsg_split_tfsm.derive_all_guard_tran_mutants(split_tran1)
                        tran_mutants += tran_guard_mutants
                        tran_guard_mutants = tsg_split_tfsm.derive_all_guard_tran_mutants(split_tran2)
                        tran_mutants += tran_guard_mutants
            fo_mutants += tran_mutants
        return fo_mutants

    def derive_all_output_tran_mutants(self, spec_tran):
        output_mutants = list()
        for o_num in range(0, self.tfsm_spec.outputs_number):
            if o_num != spec_tran.output:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.output = int(o_num)
                output_mutants.append(tfsm_mut)
        return output_mutants

    def derive_all_transition_tran_mutants(self, spec_tran):
        tran_mutants = list()
        for mut_s in range(0, self.tfsm_spec.states_number):
            if mut_s != spec_tran.end_state:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.end_state = int(mut_s)
                tran_mutants.append(tfsm_mut)
        return tran_mutants

    def derive_all_delay_tran_mutants(self, spec_tran):
        delay_mutants = list()
        for mut_d in range(1, self.tfsm_spec.d):
            if mut_d != spec_tran.delay:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.delay = int(mut_d)
                delay_mutants.append(tfsm_mut)
        return delay_mutants

    def derive_all_guard_tran_mutants(self, spec_tran):
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
        fo_tests = list()
        for tran_name in self.tfsm_spec.tfsm.transition_dict.keys():
            #print("tran_name =", tran_name)
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
                    ######
                    tsg_split_tfsm = TestSuiteGeneration(spec_split_tfsm)
                    tran_guard_mutants, tran_guard_tests = tsg_split_tfsm.derive_guard_tran_mutants_for_time_evaluation(split_tran1)
                    fo_mutants += tran_guard_mutants
                    tran_tests += tran_guard_tests
                    tran_guard_mutants, tran_guard_tests = tsg_split_tfsm.derive_guard_tran_mutants_for_time_evaluation(split_tran2)
                    fo_mutants += tran_guard_mutants
                    tran_tests += tran_guard_tests
                    self.fo_mutants_number += int(tsg_split_tfsm.fo_mutants_number)
                    '''
                    if spec_split_tfsm.tfsm.is_race_free():
                        tsg_split_tfsm = TestSuiteGeneration(spec_split_tfsm)
                        tran_guard_mutants, tran_guard_tests = tsg_split_tfsm.derive_guard_tran_mutants(split_tran1)
                        fo_mutants += tran_guard_mutants
                        tran_tests += tran_guard_tests
                        tran_guard_mutants, tran_guard_tests = tsg_split_tfsm.derive_guard_tran_mutants(split_tran2)
                        fo_mutants += tran_guard_mutants
                        tran_tests += tran_guard_tests
                    '''
            fo_mutants += tran_mutants
            fo_tests += tran_tests
        return fo_mutants, fo_tests

    def derive_guard_tran_mutants_for_time_evaluation(self, spec_tran):
        guard_mutants = list()
        guard_test_sequences = list()
        max_len = self.tfsm_spec.tfsm.ell
        for o_num in range(0, self.tfsm_spec.outputs_number):
            if o_num != spec_tran.output:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.output = int(o_num)
                untimed_seq = list(self.tfsm_spec.fsm.transfer_sequences[spec_tran.start_state])
                tis_ds = self.extend_to_timed_test_sequence(untimed_seq, "left")
                tis_ds.append_to_timed_sequence(TimedSequence([(spec_tran.input, spec_tran.time_guard[0])]))
                out_spec = self.tfsm_spec.derive_output_sequence(0, tis_ds)
                out_mut = tfsm_mut.derive_output_sequence(0, tis_ds)
                if out_spec and out_mut and out_spec != out_mut:
                    self.fo_mutants_number += 1
                    guard_mutants.append(tfsm_mut)
                    guard_test_sequences.append(tis_ds)
        for mut_s in range(0, self.tfsm_spec.states_number):
            if mut_s != spec_tran.end_state:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.end_state = int(mut_s)
                tis_ds = self.extend_to_timed_test_sequence(self.tfsm_spec.fsm.transfer_sequences[spec_tran.start_state], "left")
                tis_ds.append_to_timed_sequence(TimedSequence([(spec_tran.input, spec_tran.time_guard[0])]))
                if spec_tran.end_state < mut_s:
                    untimed_tail = list(self.tfsm_spec.fsm.shortest_DSs[spec_tran.end_state, mut_s])
                else:
                    untimed_tail = list(self.tfsm_spec.fsm.shortest_DSs[mut_s, spec_tran.end_state])
                timed_tail = self.extend_to_timed_test_sequence(untimed_tail, "left")
                tis_ds.append_to_timed_sequence(timed_tail)
                out_spec = self.tfsm_spec.derive_output_sequence(0, tis_ds)
                out_mut = tfsm_mut.derive_output_sequence(0, tis_ds)
                if out_spec and out_mut and out_spec != out_mut:
                    self.fo_mutants_number += 1
                    guard_mutants.append(tfsm_mut)
                    guard_test_sequences.append(tis_ds)
        for mut_d in range(1, self.tfsm_spec.d):
            if mut_d != spec_tran.delay:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.delay = int(mut_d)
                rf_tree = RaceFree_TruncatedTree_for_testing(self.tfsm_spec, tfsm_mut, spec_tran.start_state, mut_tran.start_state, max_len)
                tail_ds = rf_tree.derive_bfs_tree()
                if tail_ds:
                    tis_ds = self.extend_to_timed_test_sequence(self.tfsm_spec.fsm.transfer_sequences[spec_tran.start_state], "left")
                    tis_ds.append_to_timed_sequence(tail_ds)
                    out_spec = self.tfsm_spec.derive_output_sequence(0, tis_ds)
                    out_mut = tfsm_mut.derive_output_sequence(0, tis_ds)
                    if out_spec and out_mut and out_spec != out_mut:
                        self.fo_mutants_number += 1
                        guard_mutants.append(tfsm_mut)
                        guard_test_sequences.append(tis_ds)
        return guard_mutants, guard_test_sequences

    def derive_output_tran_mutants(self, spec_tran):
        output_mutants = list()
        output_test_sequences = list()
        for o_num in range(0, self.tfsm_spec.outputs_number):
            if o_num != spec_tran.output:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.output = int(o_num)
                max_len = self.tfsm_spec.states_number + self.tfsm_spec.tfsm.ell
                rf_tree = RaceFree_TruncatedTree_for_testing(self.tfsm_spec, tfsm_mut, 0, 0, max_len)
                tis_ds = rf_tree.derive_bfs_tree()
                if tis_ds:
                    output_mutants.append(tfsm_mut)
                    output_test_sequences.append(tis_ds)
        return output_mutants, output_test_sequences

    def derive_transition_tran_mutants(self, spec_tran):
        tran_mutants = list()
        tran_test_sequences = list()
        for mut_s in range(0, self.tfsm_spec.states_number):
            if mut_s != spec_tran.end_state:
                tfsm_mut = copy.deepcopy(self.tfsm_spec)
                mut_tran = tfsm_mut.tfsm.transition_dict[spec_tran.transition_name]
                mut_tran.end_state = int(mut_s)
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
                max_len = self.tfsm_spec.states_number + self.tfsm_spec.tfsm.ell
                rf_tree = RaceFree_TruncatedTree_for_testing(self.tfsm_spec, tfsm_mut, 0, 0, max_len)
                tis_ds = rf_tree.derive_bfs_tree()
                if tis_ds:
                    delay_mutants.append(tfsm_mut)
                    delay_test_sequences.append(tis_ds)
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
                    max_len = self.tfsm_spec.states_number + self.tfsm_spec.tfsm.ell
                    rf_tree = RaceFree_TruncatedTree_for_testing(self.tfsm_spec, tfsm_mut, 0, 0, max_len)
                    tis_ds = rf_tree.derive_bfs_tree()
                    if tis_ds:
                        guard_mutants.append(tfsm_mut)
                        guard_test_sequences.append(tis_ds)
        return guard_mutants, guard_test_sequences

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
                    max_len = self.tfsm_spec.tfsm.ell
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
        max_len = self.tfsm_spec.states_number + self.tfsm_spec.tfsm.ell
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
                    max_len = self.tfsm_spec.states_number + self.tfsm_spec.tfsm.ell
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
            #print("tran_name =", tran_name)
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
                        max_len = self.tfsm_spec.states_number + self.tfsm_spec.tfsm.ell
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
                        max_len = self.tfsm_spec.states_number + self.tfsm_spec.tfsm.ell
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