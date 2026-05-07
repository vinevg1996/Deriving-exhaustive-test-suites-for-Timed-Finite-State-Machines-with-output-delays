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
import statistics
import sys
import time

class EXP:
    def __init__(self, dir_name, fsm_number, u, v, d):
        #self.dir_name = "tests/fsms/"
        print("dir_name =", dir_name)
        self.dir_name = str(dir_name)
        self.fsm_number = int(fsm_number)
        self.u = int(u)
        self.v = int(v)
        self.d = int(d)
        return

    def derive_fsm(self, k):
        fsm_file = str(self.dir_name) + str(k) + '.fsm'
        fsm = FSM(fsm_file)
        fsm.parse_fsm()
        fsm.find_reachable_seq()
        fsm.derive_reverse_fsm()
        fsm.derive_shortest_DSs_bottom_up()
        return fsm

    def derive_tfsm(self, fsm):
        tfsm = RaceFreeTFSMs(fsm, self.u, self.v, self.d)
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

    def exp_untimed_ts_vs_tfsm_mutants(self):
        mutantsNumber = 100
        for k in range(0, self.fsm_number):
            print("k =", k)
            print("PRINT_FSM")
            fsm_spec = self.derive_fsm(k)
            tfsm_spec = self.derive_tfsm(fsm_spec)
            tsg_tfsm = TestSuiteGeneration(tfsm_spec)
            delay_mutants_tfsm, delay_test_tfsm = tsg_tfsm.derive_fo_mutants("delay", mutantsNumber)
            print("delay_mutants_tfsm_number =", len(delay_mutants_tfsm))
            guard_mutants_tfsm, guard_test_tfsm = tsg_tfsm.derive_fo_mutants("guard", mutantsNumber)
            print("guard_mutants_tfsm =", len(guard_mutants_tfsm))
            #{ 0, 1, 10, 11, 100, 101 }
            #ts_llm = [[0], [1], [1,0], [1,1], [1,0,0], [1,0,1]]
            ts_llm = [[0, 0, 1, 0, 1, 0, 0], [0]]
            # mean strategy
            ts_llm_mean = tsg_tfsm.extend_to_timed_test_suite(ts_llm, "mean")
            llm_mean_vs_delay = tsg_tfsm.apply_ts_to_mutants(ts_llm_mean, delay_mutants_tfsm)
            print("llm_mean_vs_delay=", llm_mean_vs_delay)
            llm_mean_vs_guard = tsg_tfsm.apply_ts_to_mutants(ts_llm_mean, guard_mutants_tfsm)
            print("llm_mean_vs_guard=", llm_mean_vs_guard)
            # mean strategy
            ts_llm_random = tsg_tfsm.extend_to_timed_test_suite(ts_llm, "random")
            llm_random_vs_delay = tsg_tfsm.apply_ts_to_mutants(ts_llm_random, delay_mutants_tfsm)
            print("llm_random_vs_delay=", llm_random_vs_delay)
            llm_random_vs_guard = tsg_tfsm.apply_ts_to_mutants(ts_llm_random, guard_mutants_tfsm)
            print("llm_random_vs_guard=", llm_random_vs_guard)
        return

    def generate_timed_mutants(self, tsg_tfsm, mutants_type, mutants_number, order):
        non_eq_mutants_list = []
        ds_mutants_list = []
        for it in range(0, mutants_number):
            #print("it =", it)
            mut, mut_ds = tsg_tfsm.derive_rand_timed_mutant(order, mutants_type)
            if mut_ds:
                non_eq_mutants_list.append(mut)
                ds_mutants_list.append(mut_ds)
        return non_eq_mutants_list, ds_mutants_list

    def generate_timed_mutants_fast(self, tsg_tfsm, mutants_number, order):
        non_eq_mutants_list = []
        for it in range(0, mutants_number):
            #print("it =", it)
            mut = tsg_tfsm.derive_rand_timed_mutant_fast(order)
            if mut:
                non_eq_mutants_list.append(mut)
        return non_eq_mutants_list

    def fault_coverage_exp(self, mutants_number, max_order):
        for k in range(0, self.fsm_number):
            #print("k =", k)
            #print("PRINT_FSM")
            fsm_spec = self.derive_fsm(k)
            fsm_spec.parse_fsm()
            # DS_based_fsm_test
            fsm_spec.find_reachable_seq()
            fsm_spec.derive_reverse_fsm()
            fsm_spec.derive_shortest_DSs_bottom_up()
            tsg = FSM_FirstOrderTestSuiteGeneration(fsm_spec)
            output_mutants, output_tss = tsg.derive_first_order_output_mutants()
            tran_mutants, tran_tss = tsg.derive_first_order_transition_mutants()
            fo_mutants = output_mutants + tran_mutants
            fo_ts = output_tss + tran_tss
            ds_fsm_ts = tsg.derive_complete_test_gradient(fo_mutants, fo_ts)
            # TT_fsm_test
            tfsm_spec = self.derive_tfsm(fsm_spec)
            ttt = TimedTransitionTour(tfsm_spec)
            ttt.derive_ttt_template(tfsm_spec)
            transition_tour = ttt.derive_ttt_projection(tfsm_spec)
            ####
            (_, _, _, ttt_random) = ttt.derive_left_mean_right_random_ttts(tfsm_spec.tfsm, 1)
            tsg_tfsm = TestSuiteGeneration(tfsm_spec)
            ds_fsm_ts_random = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "random")
            ts_llm = [[1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1], [0], [1, 0, 0]]
            ts_llm_random = tsg_tfsm.extend_to_timed_test_suite(ts_llm, "random")
            # Evaluation
            for curr_order in range(1, max_order + 1):
                print("---------")
                print("curr_order =", curr_order)
                print("---------")
                mixt_mutants, ds_mutants = self.generate_timed_mutants(tsg_tfsm, "mixt", mutants_number, curr_order)
                llm_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ts_llm_random, mixt_mutants)
                print("llm_random_vs_mixt=", llm_random_vs_mixt)
                ttt_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_random, mixt_mutants)
                print("ttt_random_vs_mixt=", ttt_random_vs_mixt)
                ds_fsm_ts_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_random, mixt_mutants)
                print("ds_fsm_ts_random_vs_mixt=", ds_fsm_ts_random_vs_mixt)
        return

    def fault_coverage_scalability_exp(self):
        rate_complete_ts_lens_average = dict()
        rate_complete_ts_lens_average[0.2] = []
        rate_complete_ts_lens_average[0.4] = []
        rate_complete_ts_lens_average[0.6] = []
        rate_complete_ts_lens_average[0.8] = []
        rate_complete_ts_lens_average[1.0] = []
        rate_complete_ts_fault_coverage_average = dict()
        rate_complete_ts_fault_coverage_average[0.2] = []
        rate_complete_ts_fault_coverage_average[0.4] = []
        rate_complete_ts_fault_coverage_average[0.6] = []
        rate_complete_ts_fault_coverage_average[0.8] = []
        rate_complete_ts_fault_coverage_average[1.0] = []
        for k in range(0, self.fsm_number):
            print("k =", k)
            fsm_spec = self.derive_fsm(k)
            fsm_spec.parse_fsm()
            tfsm_spec = self.derive_tfsm(fsm_spec)
            tsg_tfsm = TestSuiteGeneration(tfsm_spec)
            #all_fo_output_mutants, ds_output = tsg_tfsm.derive_fo_mutants("output")
            #all_fo_tran_mutants, ds_tran = tsg_tfsm.derive_fo_mutants("transition")
            all_fo_delay_mutants, ds_delay = tsg_tfsm.derive_fo_mutants("delay")
            all_fo_guard_mutants, ds_guard = tsg_tfsm.derive_fo_mutants("guard")
            all_fo_mutants = all_fo_delay_mutants + all_fo_guard_mutants
            mutants_ids_shuffled = list(range(len(all_fo_mutants)))
            random.shuffle(mutants_ids_shuffled)
            all_ds = ds_delay + ds_guard
            tsg_complete = TestSuiteGeneration_RaceFree(tfsm_spec)
            complete_ts = tsg_complete.derive_complete_test_gradient(all_fo_mutants, all_ds)
            #print("all_fo_mutants_number =", len(all_fo_mutants))
            complete_ts_lens = sum(len(tc.sequence) for tc in complete_ts)
            #print("complete_ts_lens =", complete_ts_lens)
            rate_complete_ts_lens_average[1.0].append(complete_ts_lens)
            percents = [0.2, 0.4, 0.6, 0.8]
            step = int(len(all_fo_mutants) * 0.2)
            curr_mixt_mutants = []
            curr_ds_test = []
            already_generated = []
            #for rate in percents:
            for j in range(0, len(percents)):
                #print("------------")
                #print("rate =", percents[j])
                #N = int(len(all_fo_mutants) * rate)
                for k in range(step*j,step*(j+1)):
                    curr_mixt_mutants.append(all_fo_mutants[mutants_ids_shuffled[k]])
                    curr_ds_test.append(all_ds[mutants_ids_shuffled[k]])
                curr_ds_test = tsg_complete.derive_complete_test_gradient(curr_mixt_mutants, curr_ds_test)
                curr_ds_test_len = sum(len(tc.sequence) for tc in curr_ds_test)
                #print("rate_complete_ts_lens =", curr_ds_test_len)
                rate_complete_ts_lens_average[percents[j]].append(curr_ds_test_len)
                rate_complete_ts_vs_all_mutants = tsg_tfsm.apply_ts_to_mutants(curr_ds_test, all_fo_mutants)
                #print("rate_complete_ts_vs_all_mutants=", rate_complete_ts_vs_all_mutants)
                rate_complete_ts_fault_coverage_average[percents[j]].append(rate_complete_ts_vs_all_mutants)
        for percent in [0.2, 0.4, 0.6, 0.8]:
            print("----------------------------")
            print("percent =", percent, end=": ")
            avg_len = statistics.mean(rate_complete_ts_lens_average[percent])
            print("avg_len =", avg_len)
            print("percent =", percent, end=": ")
            avg_fc = statistics.mean(rate_complete_ts_fault_coverage_average[percent])
            print("avg_fc =", avg_fc)
            print("----------------------------")
        print("percent=1.0", end=": ")
        avg_len = statistics.mean(rate_complete_ts_lens_average[1.0])
        print("avg_len =", avg_len)
        return

    def sdn_exp(self, fsm_file, tfsm_file, mutants_number, mutants_order):
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
        ds_fsm_ts = tsg.derive_complete_test_gradient(fo_mutants, fo_ts)
        #print("ds_fsm_ts_len =", len(ds_fsm_ts))
        ####
        tfsm_spec = RaceFreeTFSMs(fsm_spec, 1, 20, 10)
        tfsm_spec.parse_tfsm_from_file(tfsm_file)
        tsg_tfsm = TestSuiteGeneration(tfsm_spec)
        mixt_mutants, ds_mutants = self.generate_timed_mutants(tsg_tfsm, "mixt", mutants_number, 1)
        print("mixt_mutants_number =", len(mixt_mutants))
        tsg_complete = TestSuiteGeneration_RaceFree(tfsm_spec)
        ds_tfsm_ts = tsg_complete.derive_complete_test_gradient(mixt_mutants, ds_mutants)
        print("ds_tfsm_ts_len =", len(ds_tfsm_ts))
        #### ts_uppaal_mut
        #ts_uppaal_mut = [[0, 1, 0], [0, 1, 1, 0, 1, 0], [0, 1], [0, 1, 1, 0], [1, 0], [0, 1, 0, 1, 0], [0, 1, 1, 1, 0]]
        ts_uppaal_mut = [[0, 0], [0, 1], [1, 1], [1, 0], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1]]
        ts_uppaal_mut_random = tsg_tfsm.extend_to_timed_test_suite(ts_uppaal_mut, "random")
        print("ts_uppaal_mut_random_number =", len(ts_uppaal_mut_random))
        #### ts_uppaal_edge_cover
        ts_uppaal_edge_cover = [[0, 0, 1, 0], [1, 1, 0, 1]]
        ts_uppaal_edge_cover_random = tsg_tfsm.extend_to_timed_test_suite(ts_uppaal_edge_cover, "random")
        print("ts_uppaal_edge_cover_random_number=", len(ts_uppaal_edge_cover_random))
        #### ts_uppaal_loc_cover
        ts_uppaal_loc_cover = [[0], [0, 1], [1], [1, 0]]
        ts_uppaal_loc_cover_random = tsg_tfsm.extend_to_timed_test_suite(ts_uppaal_loc_cover, "random")
        print("ts_uppaal_loc_cover_random_number=", len(ts_uppaal_loc_cover_random))
        ### ttt_radom
        ttt = TimedTransitionTour(tfsm_spec)
        ttt.derive_ttt_template(tfsm_spec)
        transition_tour = ttt.derive_ttt_projection(tfsm_spec)
        (ttt_left, ttt_mean, ttt_right, ttt_random) = ttt.derive_left_mean_right_random_ttts(tfsm_spec.tfsm, 1)
        print("ttt_random_number=", len(ttt_random))
        #### df_fsm_radom
        ds_fsm_ts_left = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "left")
        ds_fsm_ts_mean = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "mean")
        ds_fsm_ts_right = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "right")
        ds_fsm_ts_random = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "random")
        print("ds_fsm_ts_random_number=", len(ds_fsm_ts_random))
        llm_test = [TimedSequence([(0, 1)]), TimedSequence([(0, 20)]), TimedSequence([(1, 1)]), TimedSequence([(1, 20)]), TimedSequence([(0, 1), (0, 12)]), TimedSequence([(0, 1), (0, 21)]), TimedSequence([(0, 1), (1, 2)]), TimedSequence([(0, 1), (1, 4)]), TimedSequence([(1, 1), (0, 2)]), TimedSequence([(1, 1), (0, 8)]), TimedSequence([(1, 1), (1, 9)]), TimedSequence([(1, 1), (1, 21)]), TimedSequence([(0, 1), (1, 2), (0, 10)]), TimedSequence([(0, 1), (1, 2), (0, 22)]), TimedSequence([(0, 1), (1, 2), (1, 13)]), TimedSequence([(0, 1), (1, 2), (1, 22)]), TimedSequence([(1, 1), (0, 2), (0, 10)]), TimedSequence([(1, 1), (0, 2), (0, 22)]), TimedSequence([(1, 1), (0, 2), (1, 13)]), TimedSequence([(1, 1), (0, 2), (1, 22)]), TimedSequence([(0, 1), (0, 12), (1, 13)]), TimedSequence([(0, 1), (0, 21), (1, 22)]), TimedSequence([(1, 1), (1, 9), (0, 10)]), TimedSequence([(1, 1), (1, 21), (0, 22)])]
        ### lens
        print("ds_tfsm_ts_lens =", sum(len(tc.sequence) for tc in ds_tfsm_ts))
        print("ts_uppaal_mut_random_lens =", sum(len(tc.sequence) for tc in ts_uppaal_mut_random))
        print("ts_uppaal_edge_cover_random_lens =", sum(len(tc.sequence) for tc in ts_uppaal_edge_cover_random))
        print("ts_uppaal_loc_cover_random_lens =", sum(len(tc.sequence) for tc in ts_uppaal_loc_cover_random))
        print("ttt_random_lens =", sum(len(tc.sequence) for tc in ttt_random))
        print("ds_fsm_ts_random_lens =", sum(len(tc.sequence) for tc in ds_fsm_ts_random))
        print("llm_test_lens =", sum(len(tc.sequence) for tc in llm_test))
        ##### fault coverage
        for curr_order in range(1, mutants_order):
            print("curr_order =", curr_order)
            print("------------------------")
            ds_tfsm_ts_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_tfsm_ts, mixt_mutants)
            print("ds_tfsm_ts_vs_mixt=", ds_tfsm_ts_vs_mixt)
            ##### ttt
            ttt_left_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_left, mixt_mutants)
            print("ttt_left_vs_mixt=", ttt_left_vs_mixt)
            ttt_mean_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_mean, mixt_mutants)
            print("ttt_mean_vs_mixt=", ttt_mean_vs_mixt)
            ttt_right_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_right, mixt_mutants)
            print("ttt_right_vs_mixt=", ttt_right_vs_mixt)
            ttt_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_random, mixt_mutants)
            print("ttt_random_vs_mixt=", ttt_random_vs_mixt)
            ##### ds_fsm
            ds_fsm_ts_left_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_left, mixt_mutants)
            print("ds_fsm_ts_left_vs_mixt=", ds_fsm_ts_left_vs_mixt)
            ds_fsm_ts_mean_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_mean, mixt_mutants)
            print("ds_fsm_ts_mean_vs_mixt=", ds_fsm_ts_mean_vs_mixt)
            ds_fsm_ts_right_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_right, mixt_mutants)
            print("ds_fsm_ts_right_vs_mixt=", ds_fsm_ts_right_vs_mixt)
            ds_fsm_ts_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_random, mixt_mutants)
            print("ds_fsm_ts_random_vs_mixt=", ds_fsm_ts_random_vs_mixt)
            ##### uppaal
            uppaal_loc_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ts_uppaal_loc_cover_random, mixt_mutants)
            print("uppaal_loc_random_vs_mixt=", uppaal_loc_random_vs_mixt)
            uppaal_edge_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ts_uppaal_edge_cover_random, mixt_mutants)
            print("uppaal_edge_random_vs_mixt=", uppaal_edge_random_vs_mixt)
            uppaal_mut_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ts_uppaal_mut_random, mixt_mutants)
            print("uppaal_mut_random_vs_mixt=", uppaal_mut_random_vs_mixt)
            ##### llm
            llm_vs_mixt = tsg_tfsm.apply_ts_to_mutants(llm_test, mixt_mutants)
            print("llm_vs_mixt=", llm_vs_mixt)
            if mutants_order - curr_order > 1:
                mixt_mutants, _ = self.generate_timed_mutants(tsg_tfsm, "mixt", mutants_number, curr_order)
        return

    def ansible_exp_4states(self, fsm_file, tfsm_file, mutants_number, max_order):
        fsm_spec = FSM(fsm_file)
        fsm_spec.parse_fsm()
        fsm_spec.find_reachable_seq()
        fsm_spec.derive_reverse_fsm()
        ds_fsm_ts = [["e0","e10","e11","e1","e2","e5","e15","e16","e6","e7","e8","e3","e4","e4"]]
        ds_fsm_ts += [["e2","e9","e4"],["e0","e12","e17","e18","e13","e14","e14"],["e0","e12","e19","e20"],["e2","e6","e5","e19"]]
        ds_fsm_ts += [["e0","e14","e10","e14"],["e0","e12","e15","e19"],["e0","e12","e17","e19"],["e0","e12","e20","e19"]]
        ds_fsm_ts += [["e1","e0","e14"],["e4","e0","e14"],["e2","e7","e5","e19"],["e2","e8","e5","e19"]]
        ds_fsm_ts += [["e2","e9","e5","e19"],["e2","e9","e5","e19"],["e0","e10","e12","e19"],["e0","e6","e0","e14"]]
        tfsm_spec = RaceFreeTFSMs(fsm_spec, 225, 1000, 797)
        tfsm_spec.parse_tfsm_from_file(tfsm_file)
        tsg_tfsm = TestSuiteGeneration(tfsm_spec)
        ds_fsm_ts_left = tsg_tfsm.transition_extend_to_timed_test_suite(ds_fsm_ts, "left")
        ds_fsm_ts_mean = tsg_tfsm.transition_extend_to_timed_test_suite(ds_fsm_ts, "mean")
        ds_fsm_ts_right = tsg_tfsm.transition_extend_to_timed_test_suite(ds_fsm_ts, "right")
        ds_fsm_ts_random = tsg_tfsm.transition_extend_to_timed_test_suite(ds_fsm_ts, "random")
        mixt_mutants, ds_mutants = self.generate_timed_mutants(tsg_tfsm, "mixt", mutants_number, 1)
        tsg_complete = TestSuiteGeneration_RaceFree(tfsm_spec)
        ds_tfsm_ts = tsg_complete.derive_complete_test_gradient(mixt_mutants, ds_mutants)
        llm_test = [TimedSequence([(0, 800)]), TimedSequence([(0, 1000)]), TimedSequence([(1, 800)]), TimedSequence([(1, 1000)]), TimedSequence([(2, 800)]), TimedSequence([(2, 1000)]), TimedSequence([(3, 800)]), TimedSequence([(3, 1000)]), TimedSequence([(4, 800)]), TimedSequence([(4, 1000)]), TimedSequence([(2, 800), (0, 1600)]), TimedSequence([(2, 800), (0, 1800)]), TimedSequence([(2, 800), (1, 1600)]), TimedSequence([(2, 800), (1, 1800)]), TimedSequence([(2, 800), (2, 1600)]), TimedSequence([(2, 800), (2, 1800)]), TimedSequence([(2, 800), (3, 1600)]), TimedSequence([(2, 800), (3, 1800)]), TimedSequence([(2, 800), (4, 1600)]), TimedSequence([(2, 800), (4, 1800)]), TimedSequence([(0, 800), (0, 1600)]), TimedSequence([(0, 800), (0, 1800)]), TimedSequence([(0, 800), (1, 1600)]), TimedSequence([(0, 800), (1, 1800)]), TimedSequence([(0, 800), (2, 1600)]), TimedSequence([(0, 800), (2, 1800)]), TimedSequence([(0, 800), (3, 1600)]), TimedSequence([(0, 800), (3, 1800)]), TimedSequence([(0, 800), (4, 1600)]), TimedSequence([(0, 800), (4, 1800)]), TimedSequence([(2, 800), (0, 1600), (0, 2047)]), TimedSequence([(2, 800), (0, 1600), (0, 2600)]), TimedSequence([(2, 800), (0, 1600), (1, 2400)]), TimedSequence([(2, 800), (0, 1600), (1, 2600)]), TimedSequence([(2, 800), (0, 1600), (2, 2170)]), TimedSequence([(2, 800), (0, 1600), (2, 2600)]), TimedSequence([(2, 800), (0, 1600), (3, 2400)]), TimedSequence([(2, 800), (0, 1600), (3, 2600)]), TimedSequence([(2, 800), (0, 1600), (4, 1825)]), TimedSequence([(2, 800), (0, 1600), (4, 1865)]), TimedSequence([(2, 800), (0, 1600), (5, 1990)]), TimedSequence([(2, 800), (0, 1600), (5, 2600)])]
        uppaal_mut = [["e0", "e11", "e1", "e1"], ["e0", "e11", "e1", "e3"], ["e0", "e11", "e1", "e4"],["e0", "e11", "e3", "e1"], ["e0", "e11", "e3", "e3"], ["e0", "e11", "e3", "e4"],["e0", "e11", "e4", "e1"], ["e0", "e11", "e4", "e3"], ["e0", "e11", "e4", "e4"]]
        uppaal_mut += [["e2", "e8", "e1", "e1"], ["e2", "e8", "e1", "e3"], ["e2", "e8", "e1", "e4"],["e2", "e8", "e3", "e1"],["e2", "e8","e3","e3"],["e2","e8","e3","e4"],["e2","e8","e4","e1"],["e2","e8","e4","e3"],["e2","e8","e4","e4"]]
        uppaal_mut += [["e0", "e10", "e10"], ["e0", "e10", "e13"], ["e0", "e10", "e14"],["e0", "e13", "e10"], ["e0", "e13", "e13"], ["e0", "e13", "e14"],["e0", "e14", "e10"], ["e0", "e14", "e13"], ["e0", "e14", "e14"]]
        uppaal_mut += [["e2", "e7", "e7"], ["e2", "e7", "e6"], ["e2", "e7", "e9"], ["e2", "e6", "e7"],["e2", "e6", "e6"], ["e2", "e6", "e9"], ["e2", "e9", "e7"], ["e2", "e9", "e6"],["e2", "e9", "e9"], ["e2", "e5"]]
        uppaal_mut += [["e0", "e12", "e15", "e15"], ["e0", "e12", "e15", "e16"], ["e0", "e12", "e15", "e17"],["e0", "e12", "e15", "e18"], ["e0", "e12", "e15", "e19"], ["e0", "e12", "e15", "e20"]]
        uppaal_mut += [["e0", "e12", "e17", "e15"], ["e0", "e12", "e17", "e16"], ["e0", "e12", "e17", "e17"],["e0", "e12", "e17", "e18"], ["e0", "e12", "e17", "e19"], ["e0", "e12", "e17", "e20"]]
        uppaal_mut += [["e0", "e12", "e19", "e15"], ["e0", "e12", "e19", "e16"], ["e0", "e12", "e19", "e17"],["e0", "e12", "e19", "e18"], ["e0", "e12", "e19", "e19"], ["e0", "e12", "e19", "e20"]]
        uppaal_mut += [["e0", "e12", "e20", "e15"], ["e0", "e12", "e20", "e16"], ["e0", "e12", "e20", "e17"],["e0", "e12", "e20", "e18"], ["e0", "e12", "e20", "e19"], ["e0", "e12", "e20", "e20"]]
        uppaal_mut += [["e0", "e12", "e16"], ["e0", "e12", "e18"]]
        uppaal_mut_random = tsg_tfsm.transition_extend_to_timed_test_suite(uppaal_mut, "random")
        uppaal_loc = [["e1", "e3", "e4", "e0", "e10", "e13", "e14", "e11", "e2", "e6", "e7", "e9", "e8"],["e2", "e5", "e16"], ["e0", "e12", "e15", "e17", "e17", "e19", "e20", "e18"]]
        uppaal_edge = [["e1", "e3", "e4", "e0", "e10", "e13", "e14", "e11", "e2", "e6", "e7", "e9", "e8"],["e2", "e5", "e16"], ["e0", "e12", "e19", "e20", "e15", "e20", "e17", "e20", "e18", ]]
        uppaal_loc_random = tsg_tfsm.transition_extend_to_timed_test_suite(uppaal_loc, "random")
        uppaal_edge_random = tsg_tfsm.transition_extend_to_timed_test_suite(uppaal_edge, "random")
        ### ttt
        ttt = TimedTransitionTour(tfsm_spec)
        ttt.derive_ttt_template(tfsm_spec)
        transition_tour = ttt.derive_ttt_projection(tfsm_spec)
        (ttt_left, ttt_mean, ttt_right, ttt_random) = ttt.derive_left_mean_right_random_ttts(tfsm_spec.tfsm, 1)
        ### lens
        print("ds_tfsm_ts_lens =", sum(len(tc.sequence) for tc in ds_tfsm_ts))
        print("ts_uppaal_mut_random_lens =", sum(len(tc.sequence) for tc in uppaal_mut_random))
        print("ts_uppaal_edge_cover_random_lens =", sum(len(tc.sequence) for tc in uppaal_edge_random))
        print("ts_uppaal_loc_cover_random_lens =", sum(len(tc.sequence) for tc in uppaal_loc_random))
        print("ttt_lens =", sum(len(tc.sequence) for tc in ttt_random))
        print("ds_fsm_ts_lens =", sum(len(tc.sequence) for tc in ds_fsm_ts_random))
        print("llm_test_lens =", sum(len(tc.sequence) for tc in llm_test))
        ##### fault coverage
        for curr_order in range(1, max_order):
            print("curr_order =", curr_order)
            print("------------------------")
            ds_tfsm_ts_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_tfsm_ts, mixt_mutants)
            print("ds_tfsm_ts_vs_mixt=", ds_tfsm_ts_vs_mixt)
            ##### ttt
            ttt_left_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_left, mixt_mutants)
            print("ttt_left_vs_mixt=", ttt_left_vs_mixt)
            ttt_mean_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_mean, mixt_mutants)
            print("ttt_mean_vs_mixt=", ttt_mean_vs_mixt)
            ttt_right_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_right, mixt_mutants)
            print("ttt_right_vs_mixt=", ttt_right_vs_mixt)
            ttt_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_random, mixt_mutants)
            print("ttt_random_vs_mixt=", ttt_random_vs_mixt)
            ##### ds_fsm
            ds_fsm_ts_left_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_left, mixt_mutants)
            print("ds_fsm_ts_left_vs_mixt=", ds_fsm_ts_left_vs_mixt)
            ds_fsm_ts_mean_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_mean, mixt_mutants)
            print("ds_fsm_ts_mean_vs_mixt=", ds_fsm_ts_mean_vs_mixt)
            ds_fsm_ts_right_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_right, mixt_mutants)
            print("ds_fsm_ts_right_vs_mixt=", ds_fsm_ts_right_vs_mixt)
            ds_fsm_ts_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_random, mixt_mutants)
            print("ds_fsm_ts_random_vs_mixt=", ds_fsm_ts_random_vs_mixt)
            ##### uppaal
            uppaal_loc_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(uppaal_loc_random, mixt_mutants)
            print("uppaal_loc_random_vs_mixt=", uppaal_loc_random_vs_mixt)
            uppaal_edge_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(uppaal_edge_random, mixt_mutants)
            print("uppaal_edge_random_vs_mixt=", uppaal_edge_random_vs_mixt)
            uppaal_mut_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(uppaal_mut_random, mixt_mutants)
            print("uppaal_mut_random_vs_mixt=", uppaal_mut_random_vs_mixt)
            ##### llm
            llm_vs_mixt = tsg_tfsm.apply_ts_to_mutants(llm_test, mixt_mutants)
            print("llm_vs_mixt=", llm_vs_mixt)
            if max_order - curr_order > 1:
                mixt_mutants, _ = self.generate_timed_mutants(tsg_tfsm, "mixt", mutants_number, curr_order)
        return

    def ansible_exp(self, fsm_file, tfsm_file, mutants_number, max_order):
        fsm_spec = FSM(fsm_file)
        fsm_spec.parse_fsm()
        fsm_spec.find_reachable_seq()
        fsm_spec.derive_reverse_fsm()
        #ds_fsm_ts = [["e0","e3","e4","e1","e2"],["e0","e5","e6"]]
        ds_fsm_ts = [["e0","e3","e6"],["e0","e4"],["e0","e5","e6"],["e0","e6","e6"],["e1"],["e2"]]
        ####
        tfsm_spec = RaceFreeTFSMs(fsm_spec, 225, 1000, 797)
        tfsm_spec.parse_tfsm_from_file(tfsm_file)
        tsg_tfsm = TestSuiteGeneration(tfsm_spec)
        ds_fsm_ts_left = tsg_tfsm.transition_extend_to_timed_test_suite(ds_fsm_ts, "left")
        ds_fsm_ts_mean = tsg_tfsm.transition_extend_to_timed_test_suite(ds_fsm_ts, "mean")
        ds_fsm_ts_right = tsg_tfsm.transition_extend_to_timed_test_suite(ds_fsm_ts, "right")
        ds_fsm_ts_random = tsg_tfsm.transition_extend_to_timed_test_suite(ds_fsm_ts, "random")
        mixt_mutants, ds_mutants = self.generate_timed_mutants(tsg_tfsm, "mixt", mutants_number, 1)
        tsg_complete = TestSuiteGeneration_RaceFree(tfsm_spec)
        ds_tfsm_ts = tsg_complete.derive_complete_test_gradient(mixt_mutants, ds_mutants)
        ### ttt
        ttt = TimedTransitionTour(tfsm_spec)
        ttt.derive_ttt_template(tfsm_spec)
        transition_tour = ttt.derive_ttt_projection(tfsm_spec)
        (ttt_left, ttt_mean, ttt_right, ttt_random) = ttt.derive_left_mean_right_random_ttts(tfsm_spec.tfsm, 1)
        uppaal_loc = [["e0", "e3", "e4", "e2"], ["e0", "e5"], ["e0", "e6"]]
        uppaal_loc_random = tsg_tfsm.transition_extend_to_timed_test_suite(uppaal_loc, "random")
        uppaal_edge = [["e0","e3","e4","e1","e2"],["e0","e5"],["e0","e6"]]
        uppaal_edge_random = tsg_tfsm.transition_extend_to_timed_test_suite(uppaal_edge, "random")
        uppaal_mut = [["e1","e1"],["e1","e2"],["e2","e1"],["e2","e2"]]
        uppaal_mut += [["e0","e3","e3"],["e0","e3","e4"],["e0","e3","e5"],["e0","e3","e6"]]
        uppaal_mut += [["e0","e4","e3"],["e0","e4","e4"],["e0","e4","e5"],["e0","e4","e6"]]
        uppaal_mut += [["e0", "e5", "e3"], ["e0", "e5", "e4"], ["e0", "e5", "e5"], ["e0", "e5", "e6"]]
        uppaal_mut += [["e0", "e6", "e3"], ["e0", "e6", "e4"], ["e0", "e6", "e5"], ["e0", "e6", "e6"]]
        uppaal_mut_random = tsg_tfsm.transition_extend_to_timed_test_suite(uppaal_mut, "random")
        llm_test = [TimedSequence([(0, 800), (0, 1247), (2, 2047)]), TimedSequence([(0, 800), (0, 1800), (2, 2600)]), TimedSequence([(0, 800), (1, 1600), (2, 2400)]), TimedSequence([(0, 800), (1, 1800), (2, 2600)]), TimedSequence([(0, 800), (2, 1025), (2, 1825)]), TimedSequence([(0, 800), (2, 1065), (2, 1865)]), TimedSequence([(0, 800), (2, 1067), (2, 1867)]), TimedSequence([(0, 800), (2, 1800), (2, 2600)]), TimedSequence([(0, 1000), (0, 1447), (2, 2247)]), TimedSequence([(0, 1000), (1, 1800), (2, 2600)]), TimedSequence([(0, 1000), (2, 1225), (2, 2025)]), TimedSequence([(0, 1000), (2, 1265), (2, 2065)]), TimedSequence([(0, 1000), (2, 1267), (2, 2067)]), TimedSequence([(0, 1000), (2, 2000), (2, 2800)]), TimedSequence([(1, 800), (0, 1600), (2, 2400)]), TimedSequence([(1, 800), (1, 1600), (2, 2400)]), TimedSequence([(1, 800), (2, 1600), (0, 2400)]), TimedSequence([(2, 800), (0, 1600), (2, 2400)]), TimedSequence([(2, 800), (1, 1600), (2, 2400)]), TimedSequence([(2, 800), (2, 1600), (0, 2400)])]
        ### lens
        print("ds_tfsm_ts_lens =", sum(len(tc.sequence) for tc in ds_tfsm_ts))
        print("ts_uppaal_mut_random_lens =", sum(len(tc.sequence) for tc in uppaal_mut_random))
        print("ts_uppaal_edge_cover_random_lens =", sum(len(tc.sequence) for tc in uppaal_edge_random))
        print("ts_uppaal_loc_cover_random_lens =", sum(len(tc.sequence) for tc in uppaal_loc_random))
        print("ttt_lens =", sum(len(tc.sequence) for tc in ttt_random))
        print("ds_fsm_ts_lens =", sum(len(tc.sequence) for tc in ds_fsm_ts_random))
        print("llm_test_lens =", sum(len(tc.sequence) for tc in llm_test))
        ##### fault coverage
        for curr_order in range(1, max_order):
            print("curr_order =", curr_order)
            print("------------------------")
            ds_tfsm_ts_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_tfsm_ts, mixt_mutants)
            print("ds_tfsm_ts_vs_mixt=", ds_tfsm_ts_vs_mixt)
            ##### ttt
            ttt_left_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_left, mixt_mutants)
            print("ttt_left_vs_mixt=", ttt_left_vs_mixt)
            ttt_mean_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_mean, mixt_mutants)
            print("ttt_mean_vs_mixt=", ttt_mean_vs_mixt)
            ttt_right_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_right, mixt_mutants)
            print("ttt_right_vs_mixt=", ttt_right_vs_mixt)
            ttt_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_random, mixt_mutants)
            print("ttt_random_vs_mixt=", ttt_random_vs_mixt)
            ##### ds_fsm
            ds_fsm_ts_left_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_left, mixt_mutants)
            print("ds_fsm_ts_left_vs_mixt=", ds_fsm_ts_left_vs_mixt)
            ds_fsm_ts_mean_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_mean, mixt_mutants)
            print("ds_fsm_ts_mean_vs_mixt=", ds_fsm_ts_mean_vs_mixt)
            ds_fsm_ts_right_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_right, mixt_mutants)
            print("ds_fsm_ts_right_vs_mixt=", ds_fsm_ts_right_vs_mixt)
            ds_fsm_ts_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_random, mixt_mutants)
            print("ds_fsm_ts_random_vs_mixt=", ds_fsm_ts_random_vs_mixt)
            ##### uppaal
            uppaal_loc_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(uppaal_loc_random, mixt_mutants)
            print("uppaal_loc_random_vs_mixt=", uppaal_loc_random_vs_mixt)
            uppaal_edge_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(uppaal_edge_random, mixt_mutants)
            print("uppaal_edge_random_vs_mixt=", uppaal_edge_random_vs_mixt)
            uppaal_mut_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(uppaal_mut_random, mixt_mutants)
            print("uppaal_mut_random_vs_mixt=", uppaal_mut_random_vs_mixt)
            ##### llm
            llm_vs_mixt = tsg_tfsm.apply_ts_to_mutants(llm_test, mixt_mutants)
            print("llm_vs_mixt=", llm_vs_mixt)
            if max_order - curr_order > 1:
                mixt_mutants, _ = self.generate_timed_mutants(tsg_tfsm, "mixt", mutants_number, curr_order)
        return

    def fault_coverage_and_length_exp(self, mutants_number, max_order):
        ds_tfsm_ts_fc = dict()
        ds_tfsm_ts_len = []
        ds_fsm_ts_fc = dict()
        ds_fsm_ts_fc["mean"] = dict()
        ds_fsm_ts_fc["random"] = dict()
        ds_fsm_ts_len = []
        ttt_ts_fc = dict()
        ttt_ts_fc["mean"] = dict()
        ttt_ts_fc["random"] = dict()
        ttt_ts_len = []
        for curr_order in range(1, max_order + 1):
            ds_tfsm_ts_fc[curr_order] = []
            ds_fsm_ts_fc["mean"][curr_order] = []
            ds_fsm_ts_fc["random"][curr_order] = []
            ttt_ts_fc["mean"][curr_order] = []
            ttt_ts_fc["random"][curr_order] = []
        for k in range(0, self.fsm_number):
            print("k =", k)
            fsm_spec = self.derive_fsm(k)
            fsm_spec.parse_fsm()
            # DS_based_fsm_test
            fsm_spec.find_reachable_seq()
            fsm_spec.derive_reverse_fsm()
            fsm_spec.derive_shortest_DSs_bottom_up()
            tsg = FSM_FirstOrderTestSuiteGeneration(fsm_spec)
            output_mutants, output_tss = tsg.derive_first_order_output_mutants()
            tran_mutants, tran_tss = tsg.derive_first_order_transition_mutants()
            fo_mutants = output_mutants + tran_mutants
            fo_ts = output_tss + tran_tss
            ds_fsm_ts = tsg.derive_complete_test_gradient(fo_mutants, fo_ts)
            # DS_based_tfsm_test
            tfsm_spec = self.derive_tfsm(fsm_spec)
            tsg_tfsm = TestSuiteGeneration(tfsm_spec)
            ds_fsm_ts_random = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "random")
            ds_fsm_ts_len.append(sum(len(tc.sequence) for tc in ds_fsm_ts_random))
            ds_fsm_ts_mean = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "mean")
            mixt_mutants, ds_mutants = self.generate_timed_mutants(tsg_tfsm, "mixt", mutants_number, 1)
            print("mixt_non_eq_mutants_number =", len(mixt_mutants))
            tsg_complete = TestSuiteGeneration_RaceFree(tfsm_spec)
            ds_tfsm_ts = tsg_complete.derive_complete_test_gradient(mixt_mutants, ds_mutants)
            ds_tfsm_ts_len.append(sum(len(tc.sequence) for tc in ds_tfsm_ts))
            # TTT
            ttt = TimedTransitionTour(tfsm_spec)
            ttt.derive_ttt_template(tfsm_spec)
            transition_tour = ttt.derive_ttt_projection(tfsm_spec)
            (_, ttt_mean, _, ttt_random) = ttt.derive_left_mean_right_random_ttts(tfsm_spec.tfsm, 1)
            ttt_ts_len.append(sum(len(tc.sequence) for tc in ttt_random))
            for curr_order in range(1, max_order+1):
                ### TFSM DS
                ds_tfsm_ts_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_tfsm_ts, mixt_mutants)
                ds_tfsm_ts_fc[curr_order].append(ds_tfsm_ts_vs_mixt)
                ### TTT
                ttt_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_random, mixt_mutants)
                ttt_ts_fc["random"][curr_order].append(ttt_random_vs_mixt)
                ttt_mean_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_mean, mixt_mutants)
                ttt_ts_fc["mean"][curr_order].append(ttt_mean_vs_mixt)
                #### FSM DS
                ds_fsm_ts_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_random, mixt_mutants)
                ds_fsm_ts_fc["random"][curr_order].append(ds_fsm_ts_random_vs_mixt)
                ds_fsm_ts_mean_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_mean, mixt_mutants)
                ds_fsm_ts_fc["mean"][curr_order].append(ds_fsm_ts_mean_vs_mixt)
                if max_order - curr_order >= 1:
                    mixt_mutants, _ = self.generate_timed_mutants(tsg_tfsm, "mixt", mutants_number, curr_order)
                    print("mixt_non_eq_mutants_number =", len(mixt_mutants))
        print("----------------")
        print("ds_tfsm_ts_len =", statistics.mean(ds_tfsm_ts_len))
        print("ds_fsm_ts_len =", statistics.mean(ds_fsm_ts_len))
        print("ttt_ts_len =", statistics.mean(ttt_ts_len))
        print("----------------")
        for curr_order in range(1, max_order + 1):
            print("curr_order =", curr_order)
            print("----------------")
            print("ds_tfsm_ts_fc =", statistics.mean(ds_tfsm_ts_fc[curr_order]))
            print("ds_fsm_mean_ts_fc =", statistics.mean(ds_fsm_ts_fc["mean"][curr_order]))
            print("ds_fsm_random_ts_fc =", statistics.mean(ds_fsm_ts_fc["random"][curr_order]))
            print("ttt_mean_ts_fc =", statistics.mean(ttt_ts_fc["mean"][curr_order]))
            print("ttt_random_ts_fc =", statistics.mean(ttt_ts_fc["random"][curr_order]))
            print("----------------")
        return

    def derive_dss(self, tfsm_spec, mut_list, max_len):
        dss_list = []
        for tfsm_mut in mut_list:
            #max_len = tfsm_spec.states_number + tfsm_spec.tfsm.ell
            rf_tree = RaceFree_TruncatedTree_for_testing(tfsm_spec, tfsm_mut, 0, 0, max_len)
            tis_ds = rf_tree.derive_bfs_tree()
            if tis_ds:
                dss_list.append(tis_ds)
            else:
                dss_list.append(TimedSequence([]))
        return dss_list

    def ds_distribution_length(self):
        DsTS_tfsm_len = []
        eq_mutants_number = []
        non_eq_mutants_number = []
        all_mutants_number = []
        ds_length_average = []
        for k in range(0, self.fsm_number):
            ds_length_distribution = dict()
            print("k =", k)
            fsm_spec = self.derive_fsm(k)
            fsm_spec.parse_fsm()
            tfsm_spec = self.derive_tfsm(fsm_spec)
            #max_len = tfsm_spec.states_number + tfsm_spec.tfsm.ell
            max_len = tfsm_spec.states_number + 1
            for ds_len in range(0, max_len + 1):
                ds_length_distribution[ds_len] = 0
            tsg_tfsm = TestSuiteGeneration(tfsm_spec)
            all_fo_output_mutants = tsg_tfsm.derive_all_fo_mutants("output")
            dss_list_output = self.derive_dss(tfsm_spec, all_fo_output_mutants, max_len)
            all_fo_tran_mutants = tsg_tfsm.derive_all_fo_mutants("transition")
            dss_list_tran = self.derive_dss(tfsm_spec, all_fo_tran_mutants, max_len)
            all_fo_delay_mutants = tsg_tfsm.derive_all_fo_mutants("delay")
            dss_list_delay = self.derive_dss(tfsm_spec, all_fo_delay_mutants, max_len)
            all_fo_guard_mutants = tsg_tfsm.derive_all_fo_mutants("guard")
            dss_list_guard = self.derive_dss(tfsm_spec, all_fo_guard_mutants, max_len)
            all_fo_mutants = all_fo_output_mutants + all_fo_tran_mutants + all_fo_delay_mutants + all_fo_guard_mutants
            dss_list = dss_list_output + dss_list_tran + dss_list_delay + dss_list_guard
            for curr_ds in dss_list:
                ds_length_distribution[len(curr_ds.sequence)] = ds_length_distribution[len(curr_ds.sequence)] + 1
            print("ds_length_distribution =", ds_length_distribution)
            all_mutants_number.append(len(dss_list))
            eq_mutants_number.append(ds_length_distribution[0])
            non_eq_mutants_number.append(len(dss_list)-ds_length_distribution[0])
            all_non_eq_mutants = []
            all_dss = []
            for curr_id in range(0, len(dss_list)):
                if len(dss_list[curr_id].sequence) > 0:
                    all_non_eq_mutants.append(all_fo_mutants[curr_id])
                    all_dss.append(dss_list[curr_id])
            tsg_complete = TestSuiteGeneration_RaceFree(tfsm_spec)
            DsTS_tfsm = tsg_complete.derive_complete_test_gradient(all_non_eq_mutants, all_dss)
            DsTS_tfsm_len.append(sum(len(tc.sequence) for tc in DsTS_tfsm))
        avg_all_mutants_number = statistics.mean(all_mutants_number)
        avg_eq_mutants_number = statistics.mean(eq_mutants_number)
        avg_non_eq_mutants_number = statistics.mean(non_eq_mutants_number)
        print("avg_non_eq_mutants_number =", avg_non_eq_mutants_number)
        percent_eq_mutants = float(avg_eq_mutants_number) / float(avg_all_mutants_number)
        print("percent_eq_mutants =", percent_eq_mutants)
        avg_DsTS_tfsm_len = statistics.mean(DsTS_tfsm_len)
        print("avg_DsTS_tfsm_len =", avg_DsTS_tfsm_len)
        return

    def untimed_fault_coverage_and_length(self, mutants_number, max_order, all_mutants_flag=True):
        ### DsTS
        ds_tfsm_ts_len = list()
        ds_tfsm_ts_fc = dict()
        ### UTS
        ds_fsm_ts_len = list()
        ds_fsm_ts_fc = dict()
        ds_fsm_ts_fc["left"] = dict()
        ds_fsm_ts_fc["mean"] = dict()
        ds_fsm_ts_fc["right"] = dict()
        ds_fsm_ts_fc["random"] = dict()
        ### TTT
        ttt_ts_len = list()
        ttt_ts_fc = dict()
        ttt_ts_fc["left"] = dict()
        ttt_ts_fc["mean"] = dict()
        ttt_ts_fc["right"] = dict()
        ttt_ts_fc["random"] = dict()
        for curr_order in range(1, max_order + 1):
            ds_tfsm_ts_fc[curr_order] = list()
            ds_fsm_ts_fc["left"][curr_order] = list()
            ds_fsm_ts_fc["mean"][curr_order] = list()
            ds_fsm_ts_fc["right"][curr_order] = list()
            ds_fsm_ts_fc["random"][curr_order] = list()
            ttt_ts_fc["left"][curr_order] = list()
            ttt_ts_fc["mean"][curr_order] = list()
            ttt_ts_fc["right"][curr_order] = list()
            ttt_ts_fc["random"][curr_order] = list()
        ### Fault coverage
        for k in range(0, self.fsm_number):
            print("k =", k)
            fsm_spec = self.derive_fsm(k)
            fsm_spec.parse_fsm()
            # DS_based_fsm_test
            fsm_spec.find_reachable_seq()
            fsm_spec.derive_reverse_fsm()
            fsm_spec.derive_shortest_DSs_bottom_up()
            tsg = FSM_FirstOrderTestSuiteGeneration(fsm_spec)
            output_mutants, output_tss = tsg.derive_first_order_output_mutants()
            tran_mutants, tran_tss = tsg.derive_first_order_transition_mutants()
            fo_mutants = output_mutants + tran_mutants
            fo_ts = output_tss + tran_tss
            ds_fsm_ts = tsg.derive_complete_test_gradient(fo_mutants, fo_ts)
            # DS_based_tfsm_test
            tfsm_spec = self.derive_tfsm(fsm_spec)
            tsg_tfsm = TestSuiteGeneration(tfsm_spec)
            ds_fsm_ts_left = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "left")
            ds_fsm_ts_mean = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "mean")
            ds_fsm_ts_right = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "right")
            ds_fsm_ts_random = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "random")
            ds_fsm_ts_len.append(sum(len(tc.sequence) for tc in ds_fsm_ts_random))
            if all_mutants_flag:
                all_fo_delay_mutants, ds_delay = tsg_tfsm.derive_fo_mutants("delay")
                all_fo_guard_mutants, ds_guard = tsg_tfsm.derive_fo_mutants("guard")
                mixt_mutants = all_fo_delay_mutants + all_fo_guard_mutants
                ds_mutants = ds_delay + ds_guard
            else:
                mixt_mutants, ds_mutants = self.generate_timed_mutants(tsg_tfsm, "mixt", mutants_number, 1)
            print("mixt_non_eq_mutants_number =", len(mixt_mutants))
            tsg_complete = TestSuiteGeneration_RaceFree(tfsm_spec)
            ds_tfsm_ts = tsg_complete.derive_complete_test_gradient(mixt_mutants, ds_mutants)
            ds_tfsm_ts_len.append(sum(len(tc.sequence) for tc in ds_tfsm_ts))
            # TTT
            ttt = TimedTransitionTour(tfsm_spec)
            ttt.derive_ttt_template(tfsm_spec)
            ttt.derive_ttt_projection(tfsm_spec)
            (ttt_left, ttt_mean, ttt_right, ttt_random) = ttt.derive_left_mean_right_random_ttts(tfsm_spec.tfsm, 1)
            ttt_ts_len.append(sum(len(tc.sequence) for tc in ttt_random))
            for curr_order in range(1, max_order+1):
                ### TFSM DS
                ds_tfsm_ts_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_tfsm_ts, mixt_mutants)
                ds_tfsm_ts_fc[curr_order].append(ds_tfsm_ts_vs_mixt)
                ### TTT
                ttt_left_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_left, mixt_mutants)
                ttt_ts_fc["left"][curr_order].append(ttt_left_vs_mixt)
                ttt_mean_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_mean, mixt_mutants)
                ttt_ts_fc["mean"][curr_order].append(ttt_mean_vs_mixt)
                ttt_right_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_right, mixt_mutants)
                ttt_ts_fc["right"][curr_order].append(ttt_right_vs_mixt)
                ttt_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_random, mixt_mutants)
                ttt_ts_fc["random"][curr_order].append(ttt_random_vs_mixt)
                #### FSM DS
                ds_fsm_ts_left_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_left, mixt_mutants)
                ds_fsm_ts_fc["left"][curr_order].append(ds_fsm_ts_left_vs_mixt)
                ds_fsm_ts_mean_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_mean, mixt_mutants)
                ds_fsm_ts_fc["mean"][curr_order].append(ds_fsm_ts_mean_vs_mixt)
                ds_fsm_ts_right_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_right, mixt_mutants)
                ds_fsm_ts_fc["right"][curr_order].append(ds_fsm_ts_right_vs_mixt)
                ds_fsm_ts_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_random, mixt_mutants)
                ds_fsm_ts_fc["random"][curr_order].append(ds_fsm_ts_random_vs_mixt)
                if max_order - curr_order >= 1:
                    mixt_mutants, _ = self.generate_timed_mutants(tsg_tfsm, "mixt", mutants_number, curr_order)
                    print("mixt_non_eq_mutants_number =", len(mixt_mutants))
        print("----------------")
        print("ds_tfsm_ts_len =", statistics.mean(ds_tfsm_ts_len))
        print("ds_fsm_ts_len =", statistics.mean(ds_fsm_ts_len))
        print("ttt_ts_len =", statistics.mean(ttt_ts_len))
        print("----------------")
        print("ds_tfsm_ts_fc = ", int(statistics.mean(ds_tfsm_ts_fc[1]) * 100), "% / ", int(statistics.mean(ds_tfsm_ts_fc[2]) * 100), "% / ", int(statistics.mean(ds_tfsm_ts_fc[3]) * 100), "%")
        print("ttt_left_ts_fc = ", int(statistics.mean(ttt_ts_fc["left"][1]) * 100), "% / ", int(statistics.mean(ttt_ts_fc["left"][2]) * 100), "% / ", int(statistics.mean(ttt_ts_fc["left"][3]) * 100), "%")
        print("ttt_mean_ts_fc = ", int(statistics.mean(ttt_ts_fc["mean"][1]) * 100), "% / ", int(statistics.mean(ttt_ts_fc["mean"][2]) * 100), "% / ", int(statistics.mean(ttt_ts_fc["mean"][3]) * 100), "%")
        print("ttt_right_ts_fc = ", int(statistics.mean(ttt_ts_fc["right"][1]) * 100), "% / ", int(statistics.mean(ttt_ts_fc["right"][2]) * 100), "% / ", int(statistics.mean(ttt_ts_fc["right"][3]) * 100), "%")
        print("ttt_random_ts_fc = ", int(statistics.mean(ttt_ts_fc["random"][1]) * 100), "% / ", int(statistics.mean(ttt_ts_fc["random"][2]) * 100), "% / ", int(statistics.mean(ttt_ts_fc["random"][3]) * 100), "%")
        print("ds_fsm_left_ts_fc = ", int(statistics.mean(ds_fsm_ts_fc["left"][1]) * 100), "% / ", int(statistics.mean(ds_fsm_ts_fc["left"][2]) * 100), "% / ", int(statistics.mean(ds_fsm_ts_fc["left"][3]) * 100), "%")
        print("ds_fsm_mean_ts_fc = ", int(statistics.mean(ds_fsm_ts_fc["mean"][1]) * 100), "% / ", int(statistics.mean(ds_fsm_ts_fc["mean"][2]) * 100), "% / ", int(statistics.mean(ds_fsm_ts_fc["mean"][3]) * 100), "%")
        print("ds_fsm_right_ts_fc = ", int(statistics.mean(ds_fsm_ts_fc["right"][1]) * 100), "% / ", int(statistics.mean(ds_fsm_ts_fc["right"][2]) * 100), "% / ", int(statistics.mean(ds_fsm_ts_fc["right"][3]) * 100), "%")
        print("ds_fsm_random_ts_fc = ", int(statistics.mean(ds_fsm_ts_fc["random"][1]) * 100), "% / ", int(statistics.mean(ds_fsm_ts_fc["random"][2]) * 100), "% / ", int(statistics.mean(ds_fsm_ts_fc["random"][3]) * 100), "%")
        print("----------------")
        return

    def generation_tfsm_test_time(self, mutants_number):
        ds_tfsm_ts_len = []
        for k in range(0, self.fsm_number):
            fsm_spec = self.derive_fsm(k)
            fsm_spec.parse_fsm()
            fsm_spec.find_reachable_seq()
            fsm_spec.derive_reverse_fsm()
            fsm_spec.derive_shortest_DSs_bottom_up()
            tfsm_spec = self.derive_tfsm(fsm_spec)
            tsg_tfsm = TestSuiteGeneration(tfsm_spec)
            all_fo_guard_mutants, ds_guard = tsg_tfsm.derive_fo_mutants("guard")
            print("all_fo_guard_mutants =", len(all_fo_guard_mutants))
            tsg_complete = TestSuiteGeneration_RaceFree(tfsm_spec)
            ds_tfsm_ts = tsg_complete.derive_complete_test_gradient(all_fo_guard_mutants, ds_guard)
            ds_tfsm_ts_len.append(sum(len(tc.sequence) for tc in ds_tfsm_ts))
        print("ds_tfsm_ts_len =", statistics.mean(ds_tfsm_ts_len))
        return

    def generation_fsm_test_time(self, mutants_number):
        fsm_spec = self.derive_fsm(2)
        fsm_spec.parse_fsm()
        fsm_spec.find_reachable_seq()
        fsm_spec.derive_reverse_fsm()
        fsm_spec.derive_shortest_DSs_bottom_up()
        tsg = FSM_FirstOrderTestSuiteGeneration(fsm_spec)
        output_mutants, output_tss = tsg.derive_first_order_output_mutants()
        tran_mutants, tran_tss = tsg.derive_first_order_transition_mutants()
        fo_mutants = output_mutants + tran_mutants
        fo_ts = output_tss + tran_tss
        ds_fsm_ts = tsg.derive_complete_test_gradient(fo_mutants, fo_ts)
        tfsm_spec = self.derive_tfsm(fsm_spec)
        tsg_tfsm = TestSuiteGeneration(tfsm_spec)
        ds_fsm_ts_left = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "left")
        ds_fsm_ts_mean = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "mean")
        ds_fsm_ts_right = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "right")
        ds_fsm_ts_random = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "random")
        return

    def generation_ttt_test_time(self, mutants_number):
        fsm_spec = self.derive_fsm(0)
        fsm_spec.parse_fsm()
        fsm_spec.find_reachable_seq()
        fsm_spec.derive_reverse_fsm()
        fsm_spec.derive_shortest_DSs_bottom_up()
        tfsm_spec = self.derive_tfsm(fsm_spec)
        ttt = TimedTransitionTour(tfsm_spec)
        ttt.derive_ttt_template(tfsm_spec)
        transition_tour = ttt.derive_ttt_projection(tfsm_spec)
        (ttt_left, ttt_mean, ttt_right, ttt_random) = ttt.derive_left_mean_right_random_ttts(tfsm_spec.tfsm, 1)
        return

    def fault_coverage_length_and_time_evaluation(self, mutants_number, max_order):
        fo_mutants_numbers = list()
        ### DsTS
        ds_tfsm_ts_len = list()
        ds_tfsm_ts_time = list()
        ds_tfsm_ts_fc = dict()
        ### UTS
        ds_fsm_ts_len = list()
        ds_fsm_ts_time = list()
        ds_fsm_ts_fc = dict()
        ds_fsm_ts_fc["left"] = dict()
        ds_fsm_ts_fc["mean"] = dict()
        ds_fsm_ts_fc["right"] = dict()
        ds_fsm_ts_fc["random"] = dict()
        ### TTT
        ttt_ts_len = list()
        ttt_ts_time = list()
        ttt_ts_fc = dict()
        ttt_ts_fc["left"] = dict()
        ttt_ts_fc["mean"] = dict()
        ttt_ts_fc["right"] = dict()
        ttt_ts_fc["random"] = dict()
        for curr_order in range(1, max_order + 1):
            ds_tfsm_ts_fc[curr_order] = list()
            ds_fsm_ts_fc["left"][curr_order] = list()
            ds_fsm_ts_fc["mean"][curr_order] = list()
            ds_fsm_ts_fc["right"][curr_order] = list()
            ds_fsm_ts_fc["random"][curr_order] = list()
            ttt_ts_fc["left"][curr_order] = list()
            ttt_ts_fc["mean"][curr_order] = list()
            ttt_ts_fc["right"][curr_order] = list()
            ttt_ts_fc["random"][curr_order] = list()
        ### Fault coverage
        for k in range(0, self.fsm_number):
            print("k =", k)
            ### UTS test generation
            fsm_spec = self.derive_fsm(k)
            fsm_spec.parse_fsm()
            start_time = time.time()
            fsm_spec.find_reachable_seq()
            fsm_spec.derive_reverse_fsm()
            fsm_spec.derive_shortest_DSs_bottom_up()
            tsg = FSM_FirstOrderTestSuiteGeneration(fsm_spec)
            output_mutants, output_tss = tsg.derive_first_order_output_mutants()
            tran_mutants, tran_tss = tsg.derive_first_order_transition_mutants()
            fo_mutants = output_mutants + tran_mutants
            fo_ts = output_tss + tran_tss
            ds_fsm_ts = tsg.derive_complete_test_gradient(fo_mutants, fo_ts)
            ds_fsm_ts_time.append(time.time() - start_time)
            ds_fsm_ts_len.append(sum(len(tc) for tc in ds_fsm_ts))
            tfsm_spec = self.derive_tfsm(fsm_spec)
            tsg_tfsm = TestSuiteGeneration(tfsm_spec)
            ds_fsm_ts_left = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "left")
            ds_fsm_ts_mean = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "mean")
            ds_fsm_ts_right = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "right")
            ds_fsm_ts_random = tsg_tfsm.extend_to_timed_test_suite(ds_fsm_ts, "random")
            ### TFSM test generation
            start_time = time.time()
            all_fo_guard_mutants, ds_guard = tsg_tfsm.derive_fo_mutants("guard")
            tsg_complete = TestSuiteGeneration_RaceFree(tfsm_spec)
            ds_tfsm_ts = tsg_complete.derive_complete_test_gradient(all_fo_guard_mutants, ds_guard)
            ds_tfsm_ts_time.append(time.time() - start_time)
            ds_tfsm_ts_len.append(sum(len(tc.sequence) for tc in ds_tfsm_ts))
            fo_mutants_numbers.append(int(tsg_tfsm.fo_mutants_number))
            ### TimedTransitionTour test generation
            start_time = time.time()
            ttt = TimedTransitionTour(tfsm_spec)
            ttt.derive_ttt_template(tfsm_spec)
            ttt.derive_ttt_projection(tfsm_spec)
            (ttt_left, ttt_mean, ttt_right, ttt_random) = ttt.derive_left_mean_right_random_ttts(tfsm_spec.tfsm, 1)
            ttt_ts_time.append(time.time() - start_time)
            ttt_ts_len.append(sum(len(tc.sequence) for tc in ttt_random))
            if len(all_fo_guard_mutants) > mutants_number:
                curr_order_mutants = all_fo_guard_mutants[0:mutants_number]
            else:
                curr_order_mutants = all_fo_guard_mutants
            for curr_order in range(1, max_order+1):
                print("curr_order =", curr_order)
                ### TFSM DS
                ds_tfsm_ts_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_tfsm_ts, curr_order_mutants)
                ds_tfsm_ts_fc[curr_order].append(ds_tfsm_ts_vs_mixt)
                ### TTT
                ttt_left_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_left, curr_order_mutants)
                ttt_ts_fc["left"][curr_order].append(ttt_left_vs_mixt)
                ttt_mean_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_mean, curr_order_mutants)
                ttt_ts_fc["mean"][curr_order].append(ttt_mean_vs_mixt)
                ttt_right_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_right, curr_order_mutants)
                ttt_ts_fc["right"][curr_order].append(ttt_right_vs_mixt)
                ttt_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ttt_random, curr_order_mutants)
                ttt_ts_fc["random"][curr_order].append(ttt_random_vs_mixt)
                #### FSM DS
                ds_fsm_ts_left_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_left, curr_order_mutants)
                ds_fsm_ts_fc["left"][curr_order].append(ds_fsm_ts_left_vs_mixt)
                ds_fsm_ts_mean_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_mean, curr_order_mutants)
                ds_fsm_ts_fc["mean"][curr_order].append(ds_fsm_ts_mean_vs_mixt)
                ds_fsm_ts_right_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_right, curr_order_mutants)
                ds_fsm_ts_fc["right"][curr_order].append(ds_fsm_ts_right_vs_mixt)
                ds_fsm_ts_random_vs_mixt = tsg_tfsm.apply_ts_to_mutants(ds_fsm_ts_random, curr_order_mutants)
                ds_fsm_ts_fc["random"][curr_order].append(ds_fsm_ts_random_vs_mixt)
                if max_order - curr_order >= 1:
                    curr_order_mutants = self.generate_timed_mutants_fast(tsg_tfsm, mutants_number, curr_order)
        print("----------------")
        print("fo_mutants_numbers =", statistics.mean(fo_mutants_numbers))
        print("ds_tfsm_ts_len     =", statistics.mean(ds_tfsm_ts_len))
        print("ds_fsm_ts_len      =", statistics.mean(ds_fsm_ts_len))
        print("ttt_ts_len         =", statistics.mean(ttt_ts_len))
        print("----------------")
        print("ds_tfsm_ts_time =", round(statistics.mean(ds_tfsm_ts_time), 4), end=' sec\n')
        print("ds_fsm_ts_time  =", round(statistics.mean(ds_fsm_ts_time), 4), end=' sec\n')
        print("ttt_ts_time     =", round(statistics.mean(ttt_ts_time), 4), end=' sec\n')
        print("----------------")
        if max_order == 3:
            print("DsTS_fc    = ", int(statistics.mean(ds_tfsm_ts_fc[1]) * 100),end='% / ')
            print(int(statistics.mean(ds_tfsm_ts_fc[2]) * 100), end='% / ')
            print(int(statistics.mean(ds_tfsm_ts_fc[3]) * 100), end='%\n')
            print("TTT_left   = ", int(statistics.mean(ttt_ts_fc["left"][1]) * 100), end='% / ')
            print(int(statistics.mean(ttt_ts_fc["left"][2]) * 100), end='% / ')
            print(int(statistics.mean(ttt_ts_fc["left"][3]) * 100), end='%\n')
            print("TTT_right  = ", int(statistics.mean(ttt_ts_fc["right"][1]) * 100), end='% / ')
            print(int(statistics.mean(ttt_ts_fc["right"][2]) * 100), end='% / ')
            print(int(statistics.mean(ttt_ts_fc["right"][3]) * 100), end='%\n')
            print("TTT_mean   = ", int(statistics.mean(ttt_ts_fc["mean"][1]) * 100), end='% / ')
            print(int(statistics.mean(ttt_ts_fc["mean"][2]) * 100), end='% / ')
            print(int(statistics.mean(ttt_ts_fc["mean"][3]) * 100), end='%\n')
            print("TTT_random = ", int(statistics.mean(ttt_ts_fc["random"][1]) * 100), end='% / ')
            print(int(statistics.mean(ttt_ts_fc["random"][2]) * 100), end='% / ')
            print(int(statistics.mean(ttt_ts_fc["random"][3]) * 100), end='%\n')
            print("UTS_left   = ", int(statistics.mean(ds_fsm_ts_fc["left"][1]) * 100), end='% / ')
            print(int(statistics.mean(ds_fsm_ts_fc["left"][2]) * 100), end='% / ')
            print(int(statistics.mean(ds_fsm_ts_fc["left"][3]) * 100), end='%\n')
            print("UTS_right  = ", int(statistics.mean(ds_fsm_ts_fc["right"][1]) * 100), end='% / ')
            print(int(statistics.mean(ds_fsm_ts_fc["right"][2]) * 100), end='% / ')
            print(int(statistics.mean(ds_fsm_ts_fc["right"][3]) * 100), end='%\n')
            print("UTS_mean   = ", int(statistics.mean(ds_fsm_ts_fc["mean"][1]) * 100), end='% / ')
            print(int(statistics.mean(ds_fsm_ts_fc["mean"][2]) * 100), end='% / ')
            print(int(statistics.mean(ds_fsm_ts_fc["mean"][3]) * 100), end='%\n')
            print("UTS_random = ", int(statistics.mean(ds_fsm_ts_fc["random"][1]) * 100), end='% / ')
            print(int(statistics.mean(ds_fsm_ts_fc["random"][2]) * 100), end='% / ')
            print(int(statistics.mean(ds_fsm_ts_fc["random"][3]) * 100), end='%\n')
        print("----------------")
        return

def experiments_with_random_TFSMs():
    print("type <states_number> <u> <v> <d> <mutants_number>")
    states_number = int(sys.argv[1])
    u = int(sys.argv[2])
    v = int(sys.argv[3])
    d = int(sys.argv[4])
    mutants_number = int(sys.argv[5])
    dir_name = "random_fsms/fsms"+str(states_number)+"/"
    exp = EXP(dir_name, 10, u, v, d)
    exp.fault_coverage_length_and_time_evaluation(mutants_number, 3)
    return

def read_five_ints():
    while True:
        text = input("Type <states_number> <u> <v> <d> <mutants_number>: ")
        parts = text.split()
        if len(parts) != 5:
            print("Please type exactly 5 numbers.")
            continue
        try:
            states_number, u, v, d, mutants_number = map(int, parts)
            return states_number, u, v, d, mutants_number
        except ValueError:
            print("Please type integers only.")

def run_experiments():
    print("What do you want to run ?")
    print("1 - Run experiments with random TFSMs")
    print("2 - Run experiments with Ansible for the same subnets")
    print("3 - Run experiments with Ansible for different subnets")
    print("4 - Run experiments with SDN")
    print("5 - Exit")
    choice = input("Your choice: ")
    if choice == "1":
        states_number, u, v, d, mutants_number = read_five_ints()
        dir_name = "random_fsms/fsms" + str(states_number) + "/"
        exp = EXP(dir_name, 10, u, v, d)
        exp.fault_coverage_length_and_time_evaluation(mutants_number, 3)
    elif choice == "2":
        mutants_number = int(input("Type <mutants_number>: "))
        exp = EXP("case_studies/ansible/ansible_same_subnet", 1, 225, 1000, 797)
        exp.ansible_exp("case_studies/ansible/ansible_same_subnet/ansible_fsm_2states.fsm","case_studies/ansible/ansible_same_subnet/ansible_tfsm_2states.tfsm", mutants_number, 4)
    elif choice == "3":
        mutants_number = int(input("Type <mutants_number>: "))
        exp = EXP("case_studies/ansible/ansible_different_subnets", 1, 225, 1000, 797)
        exp.ansible_exp("case_studies/ansible/ansible_different_subnets/ansible_fsm_4states.fsm","case_studies/ansible/ansible_different_subnets/ansible_tfsm_4states.tfsm", mutants_number, 4)
    elif choice == "4":
        mutants_number = int(input("Type <mutants_number>: "))
        exp = EXP("case_studies/sdn", 1, 1, 20, 10)
        exp.sdn_exp("case_studies/sdn/sdn_fsm.fsm", "case_studies/sdn/sdn_tfsm.tfsm", mutants_number, 4)
    return

run_experiments()