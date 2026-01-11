from tfsm import *

class TimedTransitionTour:
    def __init__(self, my_tfsm):
        self.template_ttt = list()
        self.covered_trans = dict()
        self.set_sequence = dict()
        self.is_reached_state = dict()
        for state in my_tfsm.tfsm.tfsm.keys():
            self.set_sequence[state] = list()
            self.is_reached_state[state] = False
        for tr_name in my_tfsm.tfsm.transition_dict.keys():
            self.covered_trans[tr_name] = False
        return

    def is_there_not_reached_state(self, my_tfsm):
        for state in self.is_reached_state.keys():
            if not(self.is_reached_state[state]):
                return True
        return False

    def return_not_reached_state(self, my_tfsm):
        for state in self.is_reached_state.keys():
            if not(self.is_reached_state[state]):
                return state
        return "error"

    def is_there_not_covered_transition(self, my_tfsm):
        for tr_name in my_tfsm.tfsm.transition_dict.keys():
            if not (self.covered_trans[tr_name]):
                return True
        return False

    def return_not_covered_transition(self, my_tfsm):
        for tr_name in my_tfsm.transition_dict.keys():
            if not (self.covered_trans[tr_name]):
                return my_tfsm.transition_dict[tr_name]
        return "error"

    def return_next_tran_for_state(self, my_tfsm, state):
        for input in my_tfsm.tfsm[state].keys():
            for timed_guard in my_tfsm.tfsm[state][input].keys():
                curr_tr = my_tfsm.tfsm[state][input][timed_guard]
                if not self.covered_trans[curr_tr.transition_name]:
                    return curr_tr
        return None

    def derive_set_seq(self, my_tfsm):
        queue = [(int(my_tfsm.initial_state), [])]
        self.is_reached_state[my_tfsm.initial_state] = True
        self.set_sequence[my_tfsm.initial_state] = []
        while len(queue) > 0:
            (state, seq) = queue.pop()
            for input in my_tfsm.tfsm.tfsm[state].keys():
                for timed_guard in my_tfsm.tfsm.tfsm[state][input].keys():
                    tran = my_tfsm.tfsm.tfsm[state][input][timed_guard]
                    if not self.is_reached_state[tran.end_state]:
                        new_seq = list(seq)
                        new_seq.append(tran.transition_name)
                        queue.append((tran.end_state, new_seq))
                        self.is_reached_state[tran.end_state] = True
                        self.set_sequence[tran.end_state] = list(new_seq)
        return

    def derive_ttt_template(self, my_tfsm):
        self.derive_set_seq(my_tfsm)
        while self.is_there_not_covered_transition(my_tfsm):
            template_seq = list()
            tran = self.return_not_covered_transition(my_tfsm.tfsm)
            for tr_name in self.set_sequence[tran.start_state]:
                template_seq.append(tr_name)
            template_seq.append(tran.transition_name)
            for tr_name in template_seq:
                self.covered_trans[tr_name] = True
            state = int(tran.end_state)
            deadlock_flag = False
            while not deadlock_flag:
                curr_tr = self.return_next_tran_for_state(my_tfsm.tfsm, state)
                if curr_tr is None:
                    deadlock_flag = True
                else:
                    template_seq.append(curr_tr.transition_name)
                    self.covered_trans[curr_tr.transition_name] = True
                    state = int(curr_tr.end_state)
            self.template_ttt.append(list(template_seq))
        return

    def derive_ttt_projection(self, my_tfsm):
        ttt_proj = list()
        for ttt_seq in self.template_ttt:
            ttt_proj_seq = list()
            for tran_name in ttt_seq:
                tran = my_tfsm.tfsm.transition_dict[tran_name]
                ttt_proj_seq.append(tran.input)
            ttt_proj.append(ttt_proj_seq)
        return ttt_proj

    def derive_left_mean_right_random_ttts(self, my_tfsm, step):
        ttt_left = list()
        ttt_mean = list()
        ttt_right = list()
        ttt_random = list()
        for ttt_seq in self.template_ttt:
            ttt_left_seq = TimedSequence([])
            ttt_mean_seq = TimedSequence([])
            ttt_right_seq = TimedSequence([])
            ttt_random_seq = TimedSequence([])
            for tran_name in ttt_seq:
                tran = my_tfsm.transition_dict[tran_name]
                timestamp_left = int(tran.time_guard[0]) + step
                ttt_left_seq.append_to_timed_sequence(TimedSequence([(tran.input, timestamp_left)]))
                timestamp_mean = (float(tran.time_guard[0]) + float(tran.time_guard[1])) / 2
                ttt_mean_seq.append_to_timed_sequence(TimedSequence([(tran.input, timestamp_mean)]))
                timestamp_right = int(tran.time_guard[1]) - step
                ttt_right_seq.append_to_timed_sequence(TimedSequence([(tran.input, timestamp_right)]))
                timestamp_random = random.uniform(tran.time_guard[0], tran.time_guard[1])
                ttt_random_seq.append_to_timed_sequence(TimedSequence([(tran.input, timestamp_random)]))
            ttt_left.append(ttt_left_seq)
            ttt_mean.append(ttt_mean_seq)
            ttt_right.append(ttt_right_seq)
            ttt_random.append(ttt_random_seq)
        return (ttt_left, ttt_mean, ttt_right, ttt_random)

    def print_ttt(self, ttt):
        for ttt_seq in ttt:
            for timed_input in ttt_seq.sequence:
                print(timed_input, end='')
            print()
        return