import math
from collections import deque
from itertools import permutations
import random
import secrets
from fsm import *

class Transition:
    def __init__(self, transition_name, start_state, input, time_guard, output, delay, end_state):
        self.transition_name = transition_name
        self.start_state = start_state
        self.input = input
        self.time_guard = time_guard
        self.output = output
        self.delay = delay
        self.end_state = end_state

    def print(self):
        print(self.transition_name, end="= (")
        print(self.start_state, end=", ")
        print(self.input, end=", ")
        print(self.time_guard, end=", ")
        print(self.output, end=", ")
        print(self.delay, end=", ")
        print(self.end_state, end=")")
        print()
        return

    def __repr__(self):
        return f"Transition({self.transition_name}, {self.start_state}, {self.input}, {self.time_guard}, {self.output}, {self.delay}, {self.end_state})"

class TimedSequence:
    def __init__(self, sequence):
        self.sequence = list(sequence)
        if len(self.sequence) == 0:
            self.time = 0.0
        else:
            self.time = float(self.sequence[len(self.sequence) - 1][1])
        return

    def append_to_timed_sequence(self, tis2):
        for (i, t) in tis2.sequence:
            self.sequence.append((int(i), float(self.time+t)))
        self.update_time()
        return

    def update_time(self):
        if len(self.sequence) == 0:
            self.time = 0.0
        else:
            (_, self.time) = self.sequence[len(self.sequence) - 1]
        return

    def __repr__(self):
        return f"TimedSequence({self.sequence})"

    def __hash__(self):
        return hash(str(self.sequence))

class RaceFreeTFSMs:
    def __init__(self, fsm, u, v, d):
        self.fsm = copy.deepcopy(fsm)
        self.u = int(u)
        self.v = int(v)
        self.d = int(d)
        return

    def generate_race_free_tfsm(self, is_partial=False):
        race_free_number = 0
        race_free_flag = False
        while not race_free_flag:
            self.generate_tfsm_from_fsm(is_partial)
            race_free_flag = self.tfsm.is_race_free()
            race_free_number += 1
        return

    def generate_reverse_transition_dict(self):
        self.reverse_transition_dict = {s: [] for s in range(0, self.states_number)}
        for tran_name in self.tfsm.transition_dict.keys():
            tran = self.tfsm.transition_dict[tran_name]
            self.reverse_transition_dict[tran.end_state].append(tran)
        return

    def give_reverse_paths(self, tran_name, length):
        all_paths = list()
        paths_curr = [[str(tran_name)]]
        paths_plus_one = list()
        for curr_len in range(2, length+1):
            paths_plus_one = list()
            for path in paths_curr:
                curr_tran_name = path[0]
                curr_tran = self.tfsm.transition_dict[curr_tran_name]
                if len(self.reverse_transition_dict[curr_tran.start_state]) == 0:
                    all_paths.append(list(path))
                for prev_tran in self.reverse_transition_dict[curr_tran.start_state]:
                    path_plus_one = [prev_tran.transition_name] + list(path)
                    paths_plus_one.append(path_plus_one)
            paths_curr = list(paths_plus_one)
        all_paths += list(paths_plus_one)
        return all_paths

    def give_direct_paths(self, tran_name, length):
        paths_curr = [[str(tran_name)]]
        paths_plus_one = list()
        for curr_len in range(2, length+1):
            paths_plus_one = list()
            for path in paths_curr:
                curr_tran_name = path[-1]
                curr_tran = self.tfsm.transition_dict[curr_tran_name]
                for i in self.tfsm.tfsm[curr_tran.end_state].keys():
                    for timed_guard in self.tfsm.tfsm[curr_tran.end_state][i].keys():
                        next_tran = self.tfsm.tfsm[curr_tran.end_state][i][timed_guard]
                        path_plus_one = list(path) + [next_tran.transition_name]
                        paths_plus_one.append(path_plus_one)
            paths_curr = list(paths_plus_one)
        return paths_plus_one

    def is_race_free_tran(self, tran):
        paths = self.give_direct_paths(tran.transition_name, self.tfsm.ell)
        paths += self.give_reverse_paths(tran.transition_name, self.tfsm.ell)
        for path in paths:
            (path_race_free_flag, i, j) = self.tfsm.is_path_race_free(path)
            if not path_race_free_flag:
                return False
        return True

    def parse_tfsm_from_file(self, tfsm_file_name):
        self.tfsm = TFSM()
        self.states_number = int(self.fsm.states_number)
        self.inputs_number = int(self.fsm.inputs_number)
        self.outputs_number = int(self.fsm.outputs_number)
        self.initial_state = int(self.fsm.initial_state)
        self.tfsm.ell = math.ceil(self.d / self.u)
        self.tfsm.reverse_transition_dict = {}
        file = open(tfsm_file_name, 'r')
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
            if len(parts) >= 8:
                transition_name, start_state_num, input_num, u_num, v_num, output_num, delay_num, end_state_num = parts[:8]
                s = int(start_state_num)
                i = int(input_num)
                timed_guard = (int(u_num), int(v_num))
                o = int(output_num)
                delay = int(delay_num)
                n_s = int(end_state_num)
                tran = Transition(transition_name, s, i, timed_guard, o, delay, n_s)
                if s not in self.tfsm.tfsm.keys():
                    self.tfsm.tfsm[s] = dict()
                if i not in self.tfsm.tfsm[s].keys():
                    self.tfsm.tfsm[s][i] = dict()
                self.tfsm.tfsm[s][i][timed_guard] = tran
                self.tfsm.transition_dict[transition_name] = tran
        return

    def generate_tfsm_from_fsm(self, is_partial=False):
        self.tfsm = TFSM()
        self.tfsm.ell = math.ceil(self.d / self.u)
        self.states_number = int(self.fsm.states_number)
        self.inputs_number = int(self.fsm.inputs_number)
        self.outputs_number = int(self.fsm.outputs_number)
        self.initial_state = int(self.fsm.initial_state)
        self.tfsm.reverse_transition_dict = {}
        j = 0
        for s in range(0, self.fsm.states_number):
            self.tfsm.tfsm[s] = dict()
            for i in range(0, self.fsm.inputs_number):
                self.tfsm.tfsm[s][i] = dict()
                [o, n_s] = self.fsm.transition_dict[s][i]
                tran_u = int(self.u)
                tran_v = int(self.v)
                tran_d = random.randint(1, self.d)
                time_guard = (tran_u, tran_v)
                transition_name = "e" + str(j)
                tran = Transition(transition_name, s, i, time_guard, o, tran_d, n_s)
                self.tfsm.tfsm[s][i][time_guard] = tran
                self.tfsm.transition_dict[transition_name] = tran
                j += 1
        return

    def derive_output_sequence(self, state, timed_input_seq):
        current_state = int(state)
        t_previous = 0
        output_set = []
        tfsm_dict = self.tfsm.tfsm
        for input, t in timed_input_seq.sequence:
            found_transition = False
            if (current_state in tfsm_dict) and (input in tfsm_dict[current_state]):
                for time_guard in tfsm_dict[current_state][input].keys():
                    t_local = t - t_previous
                    #if (t_local > time_guard[0]) and (t_local < time_guard[1]):
                    if (t_local > time_guard[0]) and (t_local <= time_guard[1]):
                        tran = tfsm_dict[current_state][input][time_guard]
                        found_transition = True
                        output_time = t + tran.delay
                        output_set.append((tran.output, output_time))
                        current_state = tran.end_state
                        t_previous = t
            if not found_transition:
                return None
                #raise ValueError(f"No valid transition found for input {input} at time {t}")
        timed_output_seq = sorted(output_set, key=lambda x: x[1])
        output_seq = [x[0] for x in timed_output_seq]
        return output_seq

class TFSM:
    def __init__(self):
        self.tfsm = {}
        self.transition_dict = {}
        self.ell = None

    def find_reachable_seq(self):
        self.transfer_sequences = dict()
        visited_states = dict()
        for s in self.tfsm:
            self.transfer_sequences[s] = TimedSequence([])
            visited_states[s] = False
        queue = [str(self.initial_state)]
        visited_states[self.initial_state] = True
        while len(queue) > 0:
            curr_state = queue.pop(0)
            for i in self.tfsm[curr_state]:
                for (u,v) in self.tfsm[curr_state][i]:
                    next_state = str(self.tfsm[curr_state][i][(u,v)].end_state)
                    if not visited_states[next_state]:
                        visited_states[next_state] = True
                        self.transfer_sequences[next_state] = TimedSequence(self.transfer_sequences[curr_state].sequence)
                        t = float(u + v) / 2
                        timed_input = TimedSequence([(i, t)])
                        self.transfer_sequences[next_state].append_to_timed_sequence(timed_input)
                        queue.append(next_state)
        return

    def give_paths_len_plus_1(self, paths):
        paths_plus_1 = []
        for path in paths:
            tran_name = path[-1]
            tran = self.transition_dict[tran_name]
            for input in self.tfsm[tran.end_state].keys():
                for timed_guard in self.tfsm[tran.end_state][input].keys():
                    path_plus_1 = list(path)
                    next_tran = self.tfsm[tran.end_state][input][timed_guard]
                    path_plus_1.append(str(next_tran.transition_name))
                    paths_plus_1.append(path_plus_1)
        return paths_plus_1

    def give_all_paths(self, state, length):
        paths = []
        paths_curr_len = []
        for input in self.tfsm[state].keys():
            for timed_guard in self.tfsm[state][input].keys():
                tran = self.tfsm[state][input][timed_guard]
                paths_curr_len.append([str(tran.transition_name)])
        #paths = paths + [paths_curr_len]
        for curr_len in range(2, length+1):
            paths_curr_len = self.give_paths_len_plus_1(paths_curr_len)
            paths = paths + paths_curr_len
        return paths

    def is_path_race_free(self, path):
        path_race_free_flag = True
        i_violation = -1
        j_violation = -1
        for i in range(0, len(path)):
            tran_i = self.transition_dict[path[i]]
            d_i = int(tran_i.delay)
            for j in range(i+1, len(path)):
                tran_j = self.transition_dict[path[j]]
                d_j = int(tran_j.delay)
                sum_u = 0
                sum_v = 0
                for k in range(i+1, j+1):
                    tran_k = self.transition_dict[path[k]]
                    (u_k, v_k) = tran_k.time_guard
                    sum_u += u_k
                    sum_v += v_k
                if (d_i - d_j > sum_u) and (d_i - d_j < sum_v):
                    path_race_free_flag = False
                    i_violation = int(i)
                    j_violation = int(j)
        return (path_race_free_flag, i_violation, j_violation)

    def is_race_free(self):
        for state in self.tfsm.keys():
            #paths = self.give_all_paths(state, self.ell+1)
            paths = self.give_all_paths(state, self.ell)
            for path in paths:
                (path_race_free_flag, i, j) = self.is_path_race_free(path)
                if not path_race_free_flag:
                    #print("path =", path)
                    #print("(i, j) =", i, j)
                    return False
        return True

    def generate_race_free_tfsm(self, fsm_file, u, v, d, is_partial=False):
        race_free_number = 0
        race_free_flag = False
        while not race_free_flag:
            file = open(fsm_file, 'r')
            self.ell = math.ceil(d / u)
            j = 0
            for line in file:
                parts = line.strip().split()
                if not parts or parts[0] in {'F', 'S', 'I', 'O', 'P'}:
                    continue
                if parts[0] == 'N0':
                    self.initial_state = parts[1]
                    continue
                if len(parts) >= 4:
                    j += 1
                    start_state, input, output, end_state = parts[:4]
                    if is_partial:
                        tran_u = random.randint(u, v-1)
                        tran_v = random.randint(tran_u + 1, v)
                    else:
                        tran_u = int(u)
                        tran_v = int(v)
                    tran_d = random.randint(1, d)
                    if start_state not in self.tfsm:
                        self.tfsm[start_state] = {}
                    if input not in self.tfsm[start_state]:
                        self.tfsm[start_state][input] = {}
                    time_guard = (tran_u, tran_v)
                    transition_name = "e" + str(j)
                    self.tfsm[start_state][input][time_guard] = Transition(transition_name, start_state, input,
                                                                           time_guard, output, int(tran_d), end_state)
                    self.transition_dict[transition_name] = Transition(transition_name, start_state, input, time_guard,
                                                                       output, int(tran_d), end_state)
            tran_u_min = random.randint(1, j)
            transition_name = "e" + str(tran_u_min)
            tran = self.transition_dict[transition_name]
            (tran_u, tran_v) = tran.time_guard
            self.transition_dict[transition_name] = Transition(transition_name, tran.start_state, tran.input, (u, tran_v),
                                                               tran.output, int(tran.delay), tran.end_state)
            tran_v_max = random.randint(1, j)
            transition_name = "e" + str(tran_v_max)
            tran = self.transition_dict[transition_name]
            (tran_u, tran_v) = tran.time_guard
            self.transition_dict[transition_name] = Transition(transition_name, tran.start_state, tran.input,(tran_u, v),
                                                               tran.output, int(tran.delay), tran.end_state)
            tran_d_max = random.randint(1, j)
            transition_name = "e" + str(tran_d_max)
            tran = self.transition_dict[transition_name]
            self.transition_dict[transition_name] = Transition(transition_name, tran.start_state, tran.input, tran.time_guard,
                                                               tran.output, int(d), tran.end_state)
            if self.is_race_free():
                race_free_flag = True
                race_free_number += 1
        print("race_free_number =", race_free_number)
        print(self.tfsm)
        return

    def parse_tfsm_file(self, file_path):
        with open(file_path, 'r') as file:
            d_max = 0
            u_min = None
            v_max = None
            for line in file:
                parts = line.strip().split()
                if not parts or parts[0] in {'F', 'S', 'I', 'P'}:
                    continue
                if parts[0] == 'O':
                    self.output_number = int(parts[1])
                    continue
                if parts[0] == 'N0':
                    self.initial_state = parts[1]
                    continue
                if len(parts) >= 7:
                    transition_name, start_state, input, time_guard, output, delay, end_state = parts[:7]
                    if int(delay) > d_max:
                        d_max = int(delay)
                    time_guard = tuple(map(int, time_guard[1:-1].split(',')))
                    if u_min is None:
                        u_min = int(time_guard[0])
                    elif int(time_guard[0]) < u_min:
                        u_min = int(time_guard[0])
                    if v_max is None:
                        v_max = int(time_guard[1])
                    elif int(time_guard[1]) < v_max:
                        v_max = int(time_guard[1])
                    # print("u =", time_guard[0])
                    if start_state not in self.tfsm:
                        self.tfsm[start_state] = {}
                    if input not in self.tfsm[start_state]:
                        self.tfsm[start_state][input] = {}
                    # self.tfsm[start_state][input][time_guard] = (output, int(delay), end_state)
                    self.tfsm[start_state][input][time_guard] = Transition(transition_name, start_state, input,
                                                                           time_guard, output, int(delay), end_state)
                    self.transition_dict[transition_name] = Transition(transition_name, start_state, input, time_guard,
                                                                       output, int(delay), end_state)
                else:
                    raise ValueError(f"Unexpected format in line: {line}")
            self.ell = math.ceil(d_max / u_min)
            self.d = int(d_max)
            self.u = int(u_min)
            return

    def compute_run(self, timed_sequence):
        if self.initial_state is None:
            raise ValueError("Initial state of TFSM is not defined.")

        current_state = self.initial_state
        t_previous = 0
        run = []

        for input, t in timed_sequence.sequence:
            found_transition = False
            if current_state in self.tfsm and input in self.tfsm[current_state]:
                for time_guard in self.tfsm[current_state][input].keys():
                    tran = self.tfsm[current_state][input][time_guard]
                    t_local = t - t_previous
                    if t_local > time_guard[0] and t_local < time_guard[1]:
                        found_transition = True
                        curr_transition = []
                        curr_transition.append(str(current_state))
                        curr_transition.append((str(tran.input), str(t)))
                        curr_transition.append((str(tran.output), str(t + tran.delay)))
                        curr_transition.append(str(tran.end_state))
                        run.append(curr_transition)
                        current_state = tran.end_state
                        t_previous = t
            if not found_transition:
                raise ValueError(f"No valid transition found for input {input} at time {t}")
        return run

    def generate_timed_output_sequences(self, timed_outputs):
        all_permutations = permutations(timed_outputs)
        valid_sequences = set()
        for perm in all_permutations:
            if all(perm[i][1] <= perm[i + 1][1] for i in range(len(perm) - 1)):
                valid_sequences.add(perm)
        return valid_sequences

    def generate_output_sequences(self, state, timed_sequence):
        current_state = str(state)
        t_previous = 0
        output_set = []

        for input, t in timed_sequence.sequence:
            found_transition = False
            if (current_state in self.tfsm) and (input in self.tfsm[current_state]):
                for time_guard in self.tfsm[current_state][input].keys():
                    t_local = t - t_previous
                    if (t_local > time_guard[0]) and (t_local < time_guard[1]):
                        tran = self.tfsm[current_state][input][time_guard]
                        found_transition = True
                        output_time = t + tran.delay
                        output_set.append((tran.output, output_time))
                        current_state = tran.end_state
                        t_previous = t
            if not found_transition:
                raise ValueError(f"No valid transition found for input {input} at time {t}")

        # Sort the output sequence by time
        output_sequences = self.generate_timed_output_sequences(output_set)
        # output_sequence.sort(key=lambda x: x[1])
        return output_sequences

    def generate_output_projections(self, state, timed_sequence):
        output_sequences = self.generate_output_sequences(state, timed_sequence)
        output_projections = list()
        for out_seq in output_sequences:
            output_projection = list()
            for (o, tau) in out_seq:
                output_projection.append(o)
            output_projections.append(output_projection)
        return output_projections