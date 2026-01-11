import re
from platform import mac_ver

import math
from fsm import *
from tfsm import *
from tfsm_test_suite_generation import *
#from tfsm_ds import *
#from tfsm_hs import *
from ttt import *
#from ds import *

def demo_fsm():
    #fsm = FSM("tests/fsms/0.fsm")
    fsm = FSM("tests/fsms/0.fsm")
    fsm.parse_fsm()
    fsm.derive_reverse_fsm()
    fsm.derive_shortest_DSs_bottom_up()
    return

def demo_FO_test_fsm():
    fsm = FSM("tests/fsms/0.fsm")
    fsm.parse_fsm()
    tsg = FSM_FirstOrderTestSuiteGeneration(fsm)
    output_mutants, output_tss = tsg.derive_first_order_output_mutants()
    tran_mutants, tran_tss = tsg.derive_first_order_transition_mutants()
    fo_mutants = output_mutants + tran_mutants
    fo_ts = output_tss + tran_tss
    complete_fsm_ts = tsg.derive_complete_test_gradient(fo_mutants, fo_ts)
    print("complete_fsm_ts =", complete_fsm_ts)
    return

def demo_tfsm():
    fsm_file = "tests/fsms/0.fsm"
    fsm = FSM(fsm_file)
    fsm.parse_fsm()
    tfsm_spec = RaceFreeTFSMs(fsm_file, 3, 4, 9)
    tfsm_spec.generate_race_free_tfsm(False)
    tis = TimedSequence([("i0", 3.5),("i0", 7.0)])
    timed_output_seq = tfsm_spec.derive_output_sequence("s0", tis)
    print("timed_output_seq =", timed_output_seq)
    return

def demo_FO_test_tfsm():
    fsm_file = "tests/fsms/0.fsm"
    fsm = FSM(fsm_file)
    fsm.parse_fsm()
    fsm.find_reachable_seq()
    fsm.derive_reverse_fsm()
    fsm.derive_shortest_DSs_bottom_up()
    # FSM FO test suites #
    tsg = FSM_FirstOrderTestSuiteGeneration(fsm)
    output_mutants, output_tss = tsg.derive_first_order_output_mutants()
    tran_mutants, tran_tss = tsg.derive_first_order_transition_mutants()
    fo_mutants = output_mutants + tran_mutants
    fo_ts = output_tss + tran_tss
    complete_fsm_ts = tsg.derive_complete_test_gradient(fo_mutants, fo_ts)
    # TFSM FO mutants #
    tfsm_spec = RaceFreeTFSMs(fsm, 3, 4, 9)
    tfsm_spec.generate_race_free_tfsm(False)
    tsg_tfsm = TestSuiteGeneration_RaceFree(tfsm_spec)
    out_muts_tfsm, out_tests_tfsm = tsg_tfsm.derive_fo_output_mutants()
    print(len(out_muts_tfsm))
    tran_muts_tfsm, tran_tests_tfsm = tsg_tfsm.derive_fo_tran_mutants()
    print(len(tran_muts_tfsm))
    delay_mutants_tfsm, delay_test_sequences_tfsm = tsg_tfsm.derive_fo_delay_mutants()
    print(len(delay_mutants_tfsm))
    fo_mutants_tfsm = out_muts_tfsm + tran_muts_tfsm + delay_mutants_tfsm
    uts = tsg_tfsm.mean_extension_for_untimed_test_suite(complete_fsm_ts)
    percent = tsg_tfsm.apply_ts_to_mutants(uts, fo_mutants_tfsm)
    print("percent =", percent)
    return

def exp_fault_coverage():
    dir = "tests/fsms/"
    fsm_number = 100
    killed_by_uts_mean_list = []
    killed_by_ttt_mean_list = []
    for k in range(0, fsm_number):
        if k % 10 == 0:
            print("k =", k)
        fsm_file = str(dir) + str(k) + '.fsm'
        fsm = FSM(fsm_file)
        fsm.parse_fsm()
        fsm.find_reachable_seq()
        fsm.derive_reverse_fsm()
        fsm.derive_shortest_DSs_bottom_up()
        tsg = FSM_FirstOrderTestSuiteGeneration(fsm)
        output_mutants, output_tss = tsg.derive_first_order_output_mutants()
        tran_mutants, tran_tss = tsg.derive_first_order_transition_mutants()
        fo_mutants = output_mutants + tran_mutants
        fo_ts = output_tss + tran_tss
        complete_fsm_ts = tsg.derive_complete_test_gradient(fo_mutants, fo_ts)
        ###############
        tfsm_spec = RaceFreeTFSMs(fsm, 3, 6, 12)
        tfsm_spec.generate_race_free_tfsm(False)
        tsg_tfsm = TestSuiteGeneration_RaceFree(tfsm_spec)
        out_muts_tfsm, out_tests_tfsm = tsg_tfsm.derive_fo_output_mutants()
        #print(len(out_muts_tfsm))
        tran_muts_tfsm, tran_tests_tfsm = tsg_tfsm.derive_fo_tran_mutants()
        #print(len(tran_muts_tfsm))
        delay_mutants_tfsm, delay_test_sequences_tfsm = tsg_tfsm.derive_fo_delay_mutants()
        #print(len(delay_mutants_tfsm))
        fo_mutants_tfsm = out_muts_tfsm + tran_muts_tfsm + delay_mutants_tfsm
        uts = tsg_tfsm.mean_extension_for_untimed_test_suite(complete_fsm_ts)
        #######
        ttt = TimedTransitionTour(tfsm_spec)
        ttt.derive_ttt_template(tfsm_spec)
        (_, ttt_mean, _, _) = ttt.derive_left_mean_right_random_ttts(tfsm_spec.tfsm, 1)
        killed_by_uts_mean_list.append(tsg_tfsm.apply_ts_to_mutants(uts, fo_mutants_tfsm))
        killed_by_ttt_mean_list.append(tsg_tfsm.apply_ts_to_mutants(ttt_mean, fo_mutants_tfsm))
    print("killed_by_uts_mean_list =", sum(killed_by_uts_mean_list) / fsm_number)
    print("killed_by_ttt_mean_list =", sum(killed_by_ttt_mean_list) / fsm_number)
    return

def derive_guard_test():
    fsm_file = "tests/fsms/0.fsm"
    fsm = FSM(fsm_file)
    fsm.parse_fsm()
    fsm.find_reachable_seq()
    fsm.derive_reverse_fsm()
    fsm.derive_shortest_DSs_bottom_up()
    tfsm_spec = RaceFreeTFSMs(fsm, 3, 5, 9)
    tfsm_spec.generate_race_free_tfsm(False)
    tsg_tfsm = TestSuiteGeneration_RaceFree(tfsm_spec)
    tsg_tfsm.derive_fo_guard_mutants()
    return

def derive_fo_mutants():
    fsm_file = "tests/fsms/0.fsm"
    fsm = FSM(fsm_file)
    fsm.parse_fsm()
    fsm.find_reachable_seq()
    fsm.derive_reverse_fsm()
    fsm.derive_shortest_DSs_bottom_up()
    tfsm_spec = RaceFreeTFSMs(fsm, 3, 5, 12)
    tfsm_spec.generate_race_free_tfsm(False)
    tsg_tfsm = TestSuiteGeneration_RaceFree(tfsm_spec)
    output_mutants, output_tss = tsg_tfsm.derive_fo_mutants("output")
    tran_mutants, tran_tss = tsg_tfsm.derive_fo_mutants("transition")
    delay_mutants, delay_tss = tsg_tfsm.derive_fo_mutants("delay")
    guard_mutants, guard_tss = tsg_tfsm.derive_fo_mutants("guard")
    print()
    return

def exp_fsm_complete_to_tfsm_complete():
    fsm_file = "tests/fsms/0.fsm"
    fsm = FSM(fsm_file)
    fsm.parse_fsm()
    fsm.find_reachable_seq()
    fsm.derive_reverse_fsm()
    fsm.derive_shortest_DSs_bottom_up()
    tsg = FSM_FirstOrderTestSuiteGeneration(fsm)
    output_mutants, output_tss = tsg.derive_first_order_output_mutants()
    tran_mutants, tran_tss = tsg.derive_first_order_transition_mutants()
    fo_mutants = output_mutants + tran_mutants
    fo_ts = output_tss + tran_tss
    complete_fsm_ts = tsg.derive_complete_test_gradient(fo_mutants, fo_ts)
    ###
    tfsm_spec = RaceFreeTFSMs(fsm, 3, 5, 12)
    tfsm_spec.generate_race_free_tfsm(False)
    return

class TestSuiteAndStatistique:
    def __init__(self, tfsm_number, test_type):
        self.test_type = str(test_type)
        self.tfsm_number = int(tfsm_number)
        self.test_suite_len = [0.0] * tfsm_number
        self.test_suite_number_of_seq = [0.0] * tfsm_number
        self.percent_of_tfsm_out_mutants_killed_by_test_suite = [0.0] * tfsm_number
        self.percent_of_tfsm_tran_mutants_killed_by_test_suite = [0.0] * tfsm_number
        self.percent_of_tfsm_delay_mutants_killed_by_test_suite = [0.0] * tfsm_number
        self.percent_of_tfsm_guard_mutants_killed_by_test_suite = [0.0] * tfsm_number
        self.percent_of_tfsm_mutants_killed_by_test_suite = [0.0] * tfsm_number
        return

    def evaluate_average(self):
        self.av_test_suite_len = float(sum(self.test_suite_len)) / self.tfsm_number
        self.av_test_suite_number_of_seq = float(sum(self.test_suite_number_of_seq)) / self.tfsm_number
        self.av_percent_of_tfsm_out_mutants_killed_by_test_suite = float(sum(self.percent_of_tfsm_out_mutants_killed_by_test_suite)) / self.tfsm_number
        self.av_percent_of_tfsm_tran_mutants_killed_by_test_suite = float(sum(self.percent_of_tfsm_tran_mutants_killed_by_test_suite)) / self.tfsm_number
        self.av_percent_of_tfsm_delay_mutants_killed_by_test_suite = float(sum(self.percent_of_tfsm_delay_mutants_killed_by_test_suite)) / self.tfsm_number
        self.av_percent_of_tfsm_guard_mutants_killed_by_test_suite = float(sum(self.percent_of_tfsm_guard_mutants_killed_by_test_suite)) / self.tfsm_number
        self.av_percent_of_tfsm_mutants_killed_by_test_suite = float(sum(self.percent_of_tfsm_mutants_killed_by_test_suite)) / self.tfsm_number
        return

    def print(self):
        print("##############")
        print(self.test_type)
        print("______________")
        print("av_test_suite_len =", self.av_test_suite_len)
        print("av_test_suite_number_of_seq =", self.av_test_suite_number_of_seq)
        print("av_percent_of_tfsm_out_mutants_killed_by_test_suite =", self.av_percent_of_tfsm_out_mutants_killed_by_test_suite)
        print("av_percent_of_tfsm_tran_mutants_killed_by_test_suite =", self.av_percent_of_tfsm_tran_mutants_killed_by_test_suite)
        print("av_percent_of_tfsm_delay_mutants_killed_by_test_suite =", self.av_percent_of_tfsm_delay_mutants_killed_by_test_suite)
        print("av_percent_of_tfsm_guard_mutants_killed_by_test_suite =", self.av_percent_of_tfsm_guard_mutants_killed_by_test_suite)
        print("av_percent_of_tfsm_mutants_killed_by_test_suite =", self.av_percent_of_tfsm_mutants_killed_by_test_suite)
        print("##############")
        return

class EXP:
    def __init__(self, fsm_number, dir):
        self.dir_name = "tests/fsms/"
        self.fsm_number = int(fsm_number)
        return

    def derive_fsm(self, k):
        fsm_file = str(self.dir_name) + str(k) + '.fsm'
        fsm = FSM(fsm_file)
        fsm.parse_fsm()
        fsm.find_reachable_seq()
        fsm.derive_reverse_fsm()
        fsm.derive_shortest_DSs_bottom_up()
        return fsm

    def derive_tfsm(self, fsm, u, v, d):
        tfsm = RaceFreeTFSMs(fsm, u, v, d)
        tfsm.generate_race_free_tfsm(False)
        return tfsm

    def derive_tt(self, tfsm):
        ttt = TimedTransitionTour(tfsm)
        ttt.derive_ttt_template(tfsm)
        transition_tour = ttt.derive_ttt_projection(tfsm)
        return transition_tour

    def derive_fd_and_test_fsm(self, fsm_spec):
        tsg = FSM_FirstOrderTestSuiteGeneration(fsm_spec)
        output_mutants, output_tss = tsg.derive_first_order_output_mutants()
        tran_mutants, tran_tss = tsg.derive_first_order_transition_mutants()
        fo_mutants = output_mutants + tran_mutants
        fo_ts = output_tss + tran_tss
        complete_fsm_ts = tsg.derive_complete_test_gradient(fo_mutants, fo_ts)
        return fo_mutants, complete_fsm_ts

    def apply_tt_to_fd(self, fsm_spec, tt, fd):
        tsg = FSM_FirstOrderTestSuiteGeneration(fsm_spec)
        percent = tsg.apply_ts_to_mutants(tt, fd)
        return percent

    def derive_fd_tfsm(self, tfsm_spec):
        tsg_tfsm = TestSuiteGeneration_RaceFree(tfsm_spec)
        out_mutants_tfsm, out_test = tsg_tfsm.derive_fo_mutants("output")
        tran_mutants_tfsm, tran_test = tsg_tfsm.derive_fo_mutants("transition")
        delay_mutants_tfsm, delay_test = tsg_tfsm.derive_fo_mutants("delay")
        guard_mutants_tfsm, guard_test = tsg_tfsm.derive_fo_mutants("guard")
        return out_mutants_tfsm, out_test, tran_mutants_tfsm, tran_test, delay_mutants_tfsm, delay_test, guard_mutants_tfsm, guard_test

    def extend_to_timed_test_suite(self, tfsm_spec, uts):
        ttt = TimedTransitionTour(tfsm_spec)
        ttt.derive_ttt_template(tfsm_spec)
        (_, ttt_mean, _, ttt_random) = ttt.derive_left_mean_right_random_ttts(tfsm_spec.tfsm, 1)
        tsg_tfsm = TestSuiteGeneration_RaceFree(tfsm_spec)
        uts_mean = tsg_tfsm.mean_extension_for_untimed_test_suite(uts)
        uts_random = tsg_tfsm.mean_extension_for_untimed_test_suite(uts, is_random=True)
        return uts_mean, uts_random, ttt_mean, ttt_random

    def apply_ts_to_timed_fd(self, tfsm_spec, out_mutants_tfsm, tran_mutants_tfsm, delay_mutants_tfsm, guard_mutants_tfsm, tts):
        tsg_tfsm = TestSuiteGeneration_RaceFree(tfsm_spec)
        fd = out_mutants_tfsm + tran_mutants_tfsm + delay_mutants_tfsm + guard_mutants_tfsm
        killed_out_mutants = tsg_tfsm.apply_ts_to_mutants(tts, out_mutants_tfsm)
        out_rate = (killed_out_mutants, len(out_mutants_tfsm))
        if len(out_mutants_tfsm):
            out_percent = float(killed_out_mutants) / len(out_mutants_tfsm)
        else:
            out_percent = 1.0
        killed_tran_mutants = tsg_tfsm.apply_ts_to_mutants(tts, tran_mutants_tfsm)
        tran_rate = (killed_tran_mutants, len(tran_mutants_tfsm))
        if len(tran_mutants_tfsm):
            tran_percent = float(killed_tran_mutants) / len(tran_mutants_tfsm)
        else:
            tran_percent = 1.0
        killed_delay_mutants = tsg_tfsm.apply_ts_to_mutants(tts, delay_mutants_tfsm)
        delay_rate = (killed_delay_mutants, len(delay_mutants_tfsm))
        if len(delay_mutants_tfsm):
            delay_percent = float(killed_delay_mutants) / len(delay_mutants_tfsm)
        else:
            delay_percent = 1.0
        killed_guard_mutants = tsg_tfsm.apply_ts_to_mutants(tts, guard_mutants_tfsm)
        guard_rate = (killed_guard_mutants, len(guard_mutants_tfsm))
        if len(guard_mutants_tfsm):
            guard_percent = float(killed_guard_mutants) / len(guard_mutants_tfsm)
        else:
            guard_percent = 1.0
        return out_percent, tran_percent, delay_percent, guard_percent, 1.0

    def exp_fsm_complete_test(self, u, v, d):
        percent_of_output_fsm_mutants_killed_by_tt = list()
        percent_of_tran_fsm_mutants_killed_by_tt = list()
        percent_of_fsm_mutants_killed_by_tt = list()
        for k in range(0, self.fsm_number):
            #print("k =", k)
            fsm_spec = self.derive_fsm(k)
            tfsm_spec = self.derive_tfsm(fsm_spec, u, v, d)
            transition_tour = self.derive_tt(tfsm_spec)
            tsg = FSM_FirstOrderTestSuiteGeneration(fsm_spec)
            output_mutants, output_tss = tsg.derive_first_order_output_mutants()
            percent_of_output_fsm_mutants_killed_by_tt.append(self.apply_tt_to_fd(fsm_spec, transition_tour, output_mutants))
            tran_mutants, tran_tss = tsg.derive_first_order_transition_mutants()
            percent_of_tran_fsm_mutants_killed_by_tt.append(self.apply_tt_to_fd(fsm_spec, transition_tour, tran_mutants))
            fd = output_mutants + tran_mutants
            percent_of_fsm_mutants_killed_by_tt.append(self.apply_tt_to_fd(fsm_spec, transition_tour, fd))
        print("------------------")
        print("states = ", fsm_spec.states_number)
        print("percent_of_output_fsm_mutants_killed_by_tt =", float(sum(percent_of_output_fsm_mutants_killed_by_tt)) / self.fsm_number)
        print("percent_of_tran_fsm_mutants_killed_by_tt =", float(sum(percent_of_tran_fsm_mutants_killed_by_tt)) / self.fsm_number)
        print("percent_of_fsm_mutants_killed_by_tt =", float(sum(percent_of_fsm_mutants_killed_by_tt)) / self.fsm_number)
        return

    def exp_fsm_complete_to_tfsm_complete_race_free(self, u, v, d):
        uts_mean_stat = TestSuiteAndStatistique(self.fsm_number, "uts_mean")
        uts_random_stat = TestSuiteAndStatistique(self.fsm_number, "uts_random")
        ttt_mean_stat = TestSuiteAndStatistique(self.fsm_number, "ttt_mean")
        ttt_random_stat = TestSuiteAndStatistique(self.fsm_number, "ttt_random")
        tfsm_complete_stat = TestSuiteAndStatistique(self.fsm_number, "tfsm_complete")
        for k in range(0, self.fsm_number):
            #if k % 10 == 0:
            #print("k =", k)
            fsm_spec = self.derive_fsm(k)
            tfsm_spec = self.derive_tfsm(fsm_spec, u, v, d)
            transition_tour = self.derive_tt(tfsm_spec)
            fo_mutants, complete_fsm_ts = self.derive_fd_and_test_fsm(fsm_spec)
            ########
            uts_mean_stat.test_suite_number_of_seq.append(len(complete_fsm_ts))
            uts_mean_stat.test_suite_len.append(sum(len(ts) for ts in complete_fsm_ts))
            ttt_mean_stat.test_suite_number_of_seq.append(len(transition_tour))
            ttt_mean_stat.test_suite_len.append(sum(len(ts) for ts in transition_tour))
            ########
            #self.percent_of_fsm_mutants_killed_by_tt.append(self.apply_tt_to_fd(fsm_spec, transition_tour, fo_mutants))
            out_mutants_tfsm, out_test_tfsm, tran_mutants_tfsm, tran_test_tfsm, delay_mutants_tfsm, delay_test_tfsm, guard_mutants_tfsm, guard_test_tfsm = self.derive_fd_tfsm(tfsm_spec)
            fd_tfsm = out_mutants_tfsm + tran_mutants_tfsm + delay_mutants_tfsm + guard_mutants_tfsm
            fo_ts_tfsm = out_test_tfsm + tran_test_tfsm + delay_test_tfsm + guard_test_tfsm
            tsg_complete = TestSuiteGeneration_RaceFree(tfsm_spec)
            complete_tfsm_ts = tsg_complete.derive_complete_test_gradient(fd_tfsm, fo_ts_tfsm)
            len_tfsm_ts = sum(len(ts.sequence) for ts in complete_tfsm_ts)
            ########
            tfsm_complete_stat.test_suite_number_of_seq.append(len(complete_tfsm_ts))
            tfsm_complete_stat.test_suite_len.append(sum(len(ts.sequence) for ts in complete_tfsm_ts))
            ########
            uts_mean, uts_random, ttt_mean, ttt_random = self.extend_to_timed_test_suite(tfsm_spec, complete_fsm_ts)
            uts_mean_vs_out, uts_mean_vs_tran, uts_mean_vs_delay, uts_mean_vs_guard, uts_mean_vs_mutants = self.apply_ts_to_timed_fd(tfsm_spec, out_mutants_tfsm, tran_mutants_tfsm, delay_mutants_tfsm, guard_mutants_tfsm, uts_mean)
            ########
            uts_mean_stat.percent_of_tfsm_out_mutants_killed_by_test_suite.append(uts_mean_vs_out)
            uts_mean_stat.percent_of_tfsm_tran_mutants_killed_by_test_suite.append(uts_mean_vs_tran)
            uts_mean_stat.percent_of_tfsm_delay_mutants_killed_by_test_suite.append(uts_mean_vs_delay)
            uts_mean_stat.percent_of_tfsm_guard_mutants_killed_by_test_suite.append(uts_mean_vs_guard)
            uts_mean_stat.percent_of_tfsm_mutants_killed_by_test_suite.append(uts_mean_vs_mutants)
            ########
            uts_random_vs_out, uts_random_vs_tran, uts_random_vs_delay, uts_random_vs_guard, uts_random_vs_mutants = self.apply_ts_to_timed_fd(tfsm_spec, out_mutants_tfsm, tran_mutants_tfsm, delay_mutants_tfsm, guard_mutants_tfsm, uts_random)
            ########
            uts_random_stat.percent_of_tfsm_out_mutants_killed_by_test_suite[k] = float(uts_random_vs_out)
            uts_random_stat.percent_of_tfsm_tran_mutants_killed_by_test_suite[k] = float(uts_random_vs_tran)
            uts_random_stat.percent_of_tfsm_delay_mutants_killed_by_test_suite[k] = float(uts_random_vs_delay)
            uts_random_stat.percent_of_tfsm_guard_mutants_killed_by_test_suite[k] = float(uts_random_vs_guard)
            uts_random_stat.percent_of_tfsm_mutants_killed_by_test_suite[k] = float(uts_random_vs_mutants)
            ########
            ttt_mean_vs_out, ttt_mean_vs_tran, ttt_mean_vs_delay, ttt_mean_vs_guard, ttt_mean_vs_mutants = self.apply_ts_to_timed_fd(tfsm_spec, out_mutants_tfsm, tran_mutants_tfsm, delay_mutants_tfsm, guard_mutants_tfsm, ttt_mean)
            ########
            ttt_mean_stat.percent_of_tfsm_out_mutants_killed_by_test_suite[k] = float(ttt_mean_vs_out)
            ttt_mean_stat.percent_of_tfsm_tran_mutants_killed_by_test_suite[k] = float(ttt_mean_vs_tran)
            ttt_mean_stat.percent_of_tfsm_delay_mutants_killed_by_test_suite[k] = float(ttt_mean_vs_delay)
            ttt_mean_stat.percent_of_tfsm_guard_mutants_killed_by_test_suite[k] = float(ttt_mean_vs_guard)
            ttt_mean_stat.percent_of_tfsm_mutants_killed_by_test_suite[k] = float(ttt_mean_vs_mutants)
            ########
            ttt_random_vs_out, ttt_random_vs_tran, ttt_random_vs_delay, ttt_random_vs_guard, ttt_random_vs_mutants = self.apply_ts_to_timed_fd(tfsm_spec, out_mutants_tfsm, tran_mutants_tfsm, delay_mutants_tfsm, guard_mutants_tfsm, ttt_random)
            ########
            ttt_random_stat.percent_of_tfsm_out_mutants_killed_by_test_suite[k] = float(ttt_random_vs_out)
            ttt_random_stat.percent_of_tfsm_tran_mutants_killed_by_test_suite[k] = float(ttt_random_vs_tran)
            ttt_random_stat.percent_of_tfsm_delay_mutants_killed_by_test_suite[k] = float(ttt_random_vs_delay)
            ttt_random_stat.percent_of_tfsm_guard_mutants_killed_by_test_suite[k] = float(ttt_random_vs_guard)
            ttt_random_stat.percent_of_tfsm_mutants_killed_by_test_suite[k] = float(ttt_random_vs_mutants)
            ########
        tfsm_complete_stat.evaluate_average()
        tfsm_complete_stat.print()
        uts_mean_stat.evaluate_average()
        uts_mean_stat.print()
        uts_random_stat.evaluate_average()
        uts_random_stat.print()
        ttt_mean_stat.evaluate_average()
        ttt_mean_stat.print()
        ttt_random_stat.evaluate_average()
        ttt_random_stat.print()
        return uts_mean_stat.av_percent_of_tfsm_mutants_killed_by_test_suite, uts_random_stat.av_percent_of_tfsm_mutants_killed_by_test_suite, ttt_mean_stat.av_percent_of_tfsm_mutants_killed_by_test_suite, ttt_random_stat.av_percent_of_tfsm_mutants_killed_by_test_suite

    def derive_fd_tfsm_all(self, tfsm_spec):
        tsg_tfsm_all = TestSuiteGeneration(tfsm_spec)
        #out_mutants_tfsm = tsg_tfsm_all.derive_fo_mutants("output")
        #tran_mutants_tfsm = tsg_tfsm_all.derive_fo_mutants("transition")
        delay_mutants_tfsm = tsg_tfsm_all.derive_fo_mutants("delay")
        guard_mutants_tfsm = tsg_tfsm_all.derive_fo_mutants("guard")
        return [], [], delay_mutants_tfsm, guard_mutants_tfsm

    def exp_fsm_complete_to_tfsm_complete_all(self, u, v, d):
        uts_mean_stat = TestSuiteAndStatistique(self.fsm_number, "uts_mean")
        uts_random_stat = TestSuiteAndStatistique(self.fsm_number, "uts_random")
        ttt_mean_stat = TestSuiteAndStatistique(self.fsm_number, "ttt_mean")
        ttt_random_stat = TestSuiteAndStatistique(self.fsm_number, "ttt_random")
        for k in range(0, self.fsm_number):
            #if k % 10 == 0:
            print("k =", k)
            fsm_spec = self.derive_fsm(k)
            tfsm_spec = self.derive_tfsm(fsm_spec, u, v, d)
            transition_tour = self.derive_tt(tfsm_spec)
            fo_mutants, complete_fsm_ts = self.derive_fd_and_test_fsm(fsm_spec)
            uts_mean_stat.test_suite_number_of_seq.append(len(complete_fsm_ts))
            uts_mean_stat.test_suite_len.append(sum(len(ts) for ts in complete_fsm_ts))
            ttt_mean_stat.test_suite_number_of_seq.append(len(transition_tour))
            ttt_mean_stat.test_suite_len.append(sum(len(ts) for ts in transition_tour))
            ########
            _, _, delay_mutants_tfsm, guard_mutants_tfsm = self.derive_fd_tfsm_all(tfsm_spec)
            uts_mean, uts_random, ttt_mean, ttt_random = self.extend_to_timed_test_suite(tfsm_spec, complete_fsm_ts)
            _, _, uts_mean_vs_delay, uts_mean_vs_guard, _ = self.apply_ts_to_timed_fd(tfsm_spec, [], [], delay_mutants_tfsm, guard_mutants_tfsm, uts_mean)
            uts_mean_stat.percent_of_tfsm_delay_mutants_killed_by_test_suite.append(uts_mean_vs_delay)
            uts_mean_stat.percent_of_tfsm_guard_mutants_killed_by_test_suite.append(uts_mean_vs_guard)
            ########
            _, _, uts_random_vs_delay, uts_random_vs_guard, _ = self.apply_ts_to_timed_fd(tfsm_spec, [], [], delay_mutants_tfsm, guard_mutants_tfsm, uts_random)
            uts_random_stat.percent_of_tfsm_delay_mutants_killed_by_test_suite[k] = float(uts_random_vs_delay)
            uts_random_stat.percent_of_tfsm_guard_mutants_killed_by_test_suite[k] = float(uts_random_vs_guard)
            ########
            _, _, ttt_mean_vs_delay, ttt_mean_vs_guard, _ = self.apply_ts_to_timed_fd(tfsm_spec, [], [], delay_mutants_tfsm, guard_mutants_tfsm, ttt_mean)
            ########
            ttt_mean_stat.percent_of_tfsm_delay_mutants_killed_by_test_suite[k] = float(ttt_mean_vs_delay)
            ttt_mean_stat.percent_of_tfsm_guard_mutants_killed_by_test_suite[k] = float(ttt_mean_vs_guard)
            _, _, ttt_random_vs_delay, ttt_random_vs_guard, _ = self.apply_ts_to_timed_fd(tfsm_spec, [], [], delay_mutants_tfsm, guard_mutants_tfsm, ttt_random)
            ttt_random_stat.percent_of_tfsm_delay_mutants_killed_by_test_suite[k] = float(ttt_random_vs_delay)
            ttt_random_stat.percent_of_tfsm_guard_mutants_killed_by_test_suite[k] = float(ttt_random_vs_guard)
        uts_mean_stat.evaluate_average()
        uts_mean_stat.print()
        uts_random_stat.evaluate_average()
        uts_random_stat.print()
        ttt_mean_stat.evaluate_average()
        ttt_mean_stat.print()
        ttt_random_stat.evaluate_average()
        ttt_random_stat.print()
        return

#demo_FO_test_fsm()
#demo_tfsm()
#demo_FO_test_tfsm()
#exp_fault_coverage()
#derive_guard_test()
#derive_fo_mutants()

def exp_with_d_increasing():
    uts_mean_fd_dict = dict()
    uts_random_dict = dict()
    ttt_mean_dict = dict()
    ttt_random_dict = dict()
    max_d = 25
    for d in range(7, max_d):
        print("---------------------------------")
        print("d =", d)
        exp = EXP(10, "tests/fsms/")
        uts_mean_fd, uts_random_fd, ttt_mean_fd, ttt_random_fd = exp.exp_fsm_complete_to_tfsm_complete_race_free(3, 5, d)
        #######
        uts_mean_fd_dict[d] = float(uts_mean_fd)
        uts_random_dict[d] = float(uts_random_fd)
        ttt_mean_dict[d] = float(ttt_mean_fd)
        ttt_random_dict[d] = float(ttt_random_fd)
        print("---------------------------------")
    print("uts_mean_fd_dict:")
    print(uts_mean_fd_dict)
    print("uts_random_dict:")
    print(uts_random_dict)
    print("ttt_mean_dict:")
    print(ttt_mean_dict)
    print("ttt_random_dict:")
    print(ttt_random_dict)
    return

def exp_with_v_minus_u_increasing():
    uts_mean_fd_dict = dict()
    uts_random_dict = dict()
    ttt_mean_dict = dict()
    ttt_random_dict = dict()
    max_d = 20
    u = 3
    max_v = max_d - 1
    for v in range(10, max_v):
        print("---------------------------------")
        print("v =", v)
        exp = EXP(10, "tests/fsms/")
        uts_mean_fd, uts_random_fd, ttt_mean_fd, ttt_random_fd = exp.exp_fsm_complete_to_tfsm_complete_race_free(u, v, max_d)
        uts_mean_fd_dict[v - u] = float(uts_mean_fd)
        uts_random_dict[v - u] = float(uts_random_fd)
        ttt_mean_dict[v - u] = float(ttt_mean_fd)
        ttt_random_dict[v - u] = float(ttt_random_fd)
        print("---------------------------------")
    print("uts_mean_fd_dict:")
    print(uts_mean_fd_dict)
    print("uts_random_dict:")
    print(uts_random_dict)
    print("ttt_mean_dict:")
    print(ttt_mean_dict)
    print("ttt_random_dict:")
    print(ttt_random_dict)
    return

def sdn_exp(fsm_file, tfsm_file):
    fsm_spec = FSM(fsm_file)
    fsm_spec.parse_fsm()
    fsm_spec.find_reachable_seq()
    fsm_spec.derive_reverse_fsm()
    fsm_spec.derive_shortest_DSs_bottom_up()
    tsg = FSM_FirstOrderTestSuiteGeneration(fsm_spec)
    output_mutants, output_tss = tsg.derive_first_order_output_mutants()
    tran_mutants, tran_tss = tsg.derive_first_order_transition_mutants()
    fo_mutants = output_mutants + tran_mutants
    fo_ts = output_tss + tran_tss
    complete_fsm_ts = tsg.derive_complete_test_gradient(fo_mutants, fo_ts)
    print("complete_fsm_ts =", complete_fsm_ts)
    ####
    tfsm_spec = RaceFreeTFSMs(fsm_spec, 1, 20, 10)
    #tfsm_spec.generate_race_free_tfsm()
    tfsm_spec.parse_tfsm_from_file(tfsm_file)
    flag = tfsm_spec.tfsm.is_race_free()
    ttt = TimedTransitionTour(tfsm_spec)
    ttt.derive_ttt_template(tfsm_spec)
    transition_tour = ttt.derive_ttt_projection(tfsm_spec)
    print("transition_tour =", transition_tour)
    ####
    (_, ttt_mean, _, ttt_random) = ttt.derive_left_mean_right_random_ttts(tfsm_spec.tfsm, 1)
    tsg_tfsm = TestSuiteGeneration_RaceFree(tfsm_spec)
    uts_mean = tsg_tfsm.mean_extension_for_untimed_test_suite(complete_fsm_ts)
    uts_random = tsg_tfsm.mean_extension_for_untimed_test_suite(complete_fsm_ts, is_random=True)
    out_mutants_tfsm, out_test_tfsm = tsg_tfsm.derive_fo_mutants("output")
    tran_mutants_tfsm, tran_test_tfsm = tsg_tfsm.derive_fo_mutants("transition")
    delay_mutants_tfsm, delay_test_tfsm = tsg_tfsm.derive_fo_mutants("delay")
    #out_mutants_tfsm, out_test_tfsm = [], []
    #tran_mutants_tfsm, tran_test_tfsm = [], []
    #delay_mutants_tfsm, delay_test_tfsm = [], []
    guard_mutants_tfsm, guard_test_tfsm = tsg_tfsm.derive_fo_mutants("guard")
    fd_tfsm = out_mutants_tfsm + tran_mutants_tfsm + delay_mutants_tfsm + guard_mutants_tfsm
    fo_ts_tfsm = out_test_tfsm + tran_test_tfsm + delay_test_tfsm + guard_test_tfsm
    tsg_complete = TestSuiteGeneration_RaceFree(tfsm_spec)
    complete_tfsm_ts = tsg_complete.derive_complete_test_gradient(fd_tfsm, fo_ts_tfsm)
    print("complete_tfsm_ts")
    for tis in complete_tfsm_ts:
        print(tis)
    return


u = 1
v = 3
d = 7
exp = EXP(10, "tests/fsms/")
exp.exp_fsm_complete_to_tfsm_complete_all(u, v, d)
#exp.exp_fsm_complete_to_tfsm_complete(u, v, d)
#sdn_exp("tests/sdn_fsm.fsm","tests/sdn_tfsm.tfsm")