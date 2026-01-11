from tfsm import *
import copy

class Node:
    def __init__(self, tfsm, s0, q0, s, q, tis, max_len, parent=None):
        self.tfsm = tfsm
        self.s0 = str(s0)
        self.q0 = str(q0)
        self.s = str(s)
        self.q = str(q)
        self.tis = copy.deepcopy(tis)
        self.max_len = int(max_len)
        self.tis_is_ds = False
        self.level = 0
        self.parent = parent
        self.children = list()
        return

    def study_status(self):
        if self.parent != None:
            self.level = int(self.parent.level) + 1
        #print("tis =", self.tis)
        out_s = self.tfsm.generate_output_projections(self.s0, self.tis)
        out_q = self.tfsm.generate_output_projections(self.q0, self.tis)
        if out_s != out_q:
            self.is_leaf = True
            self.tis_is_ds = True
        if self.level == self.max_len:
            self.is_leaf = True
        return

    def derive_successors(self, succs):
        for i in self.tfsm.tfsm[self.s].keys():
            if i in self.tfsm.tfsm[self.q].keys():
                for (u_s, v_s) in self.tfsm.tfsm[self.s][i].keys():
                    for (u_q, v_q) in self.tfsm.tfsm[self.q][i].keys():
                        tran_s = self.tfsm.tfsm[self.s][i][(u_s, v_s)]
                        next_s = str(tran_s.end_state)
                        tran_q = self.tfsm.tfsm[self.q][i][(u_q, v_q)]
                        next_q = str(tran_q.end_state)
                        u_max = max(u_s, u_q)
                        v_min = min(v_s, v_q)
                        if u_max < v_min:
                            t = float(u_max + v_min) / 2
                            time = float(self.tis.time)
                            tis_i_t = copy.deepcopy(self.tis)
                            tis_i_t.sequence.append((i, time+t))
                            tis_i_t.update_time()
                            next_node = Node(self.tfsm, self.s0, self.q0, next_s, next_q, tis_i_t, self.max_len, self)
                            #next_node = Node(self.tfsm, self.s, self.q, tis_i_t, self.max_len, self)
                            next_node.study_status()
                            self.children.append(next_node)
                            if next_node.tis_is_ds:
                                return next_node.tis
                            if next_node.level < self.max_len:
                                succs.append(next_node)
        return None

    def print_node(self):
        print("#######")
        node_str = "(" + str(self.s) + ", " + str(self.q) + ", " + str(self.level) + ")"
        print(node_str)
        print(self.tis)
        print("out_s = ", self.tfsm.generate_output_projections(self.s0, self.tis))
        print("out_q = ", self.tfsm.generate_output_projections(self.q0, self.tis))
        print("#######")
        return

class Node_for_testing:
    def __init__(self, tfsm_spec, tfsm_impl, s0, q0, s, q, tis, max_len, parent=None):
        self.tfsm_spec = tfsm_spec
        self.tfsm_impl = tfsm_impl
        self.s0 = int(s0)
        self.q0 = int(q0)
        self.s = int(s)
        self.q = int(q)
        self.tis = copy.deepcopy(tis)
        self.max_len = int(max_len)
        self.tis_is_ds = False
        self.level = 0
        self.parent = parent
        self.children = list()
        return

    def study_status(self):
        if self.parent != None:
            self.level = int(self.parent.level) + 1
        #print("tis =", self.tis)
        out_s = self.tfsm_spec.derive_output_sequence(self.s0, self.tis)
        out_q = self.tfsm_impl.derive_output_sequence(self.q0, self.tis)
        if out_s != out_q:
            self.is_leaf = True
            self.tis_is_ds = True
        if self.level == self.max_len:
            self.is_leaf = True
        return

    def derive_successors(self, succs):
        tfsm_dict_spec = self.tfsm_spec.tfsm.tfsm
        tfsm_dict_impl = self.tfsm_impl.tfsm.tfsm
        for i in tfsm_dict_spec[self.s].keys():
            for (u_s, v_s) in tfsm_dict_spec[self.s][i].keys():
                tran_s = tfsm_dict_spec[self.s][i][(u_s, v_s)]
                next_s = int(tran_s.end_state)
                if i in tfsm_dict_impl[self.q].keys():
                    for (u_q, v_q) in tfsm_dict_impl[self.q][i].keys():
                        tran_q = tfsm_dict_impl[self.q][i][(u_q, v_q)]
                        next_q = int(tran_q.end_state)
                        u_max = max(u_s, u_q)
                        v_min = min(v_s, v_q)
                        if u_max < v_min:
                            t = float(u_max + v_min) / 2
                            time = float(self.tis.time)
                            tis_i_t = copy.deepcopy(self.tis)
                            tis_i_t.sequence.append((i, time+t))
                            tis_i_t.update_time()
                            next_node = Node_for_testing(self.tfsm_spec, self.tfsm_impl, self.s0, self.q0, next_s, next_q, tis_i_t, self.max_len, self)
                            next_node.study_status()
                            self.children.append(next_node)
                            if next_node.tis_is_ds:
                                return next_node.tis
                            if next_node.level < self.max_len:
                                succs.append(next_node)
        return None

    def print_node(self):
        print("#######")
        node_str = "(" + str(self.s) + ", " + str(self.q) + ", " + str(self.level) + ")"
        print(node_str)
        print(self.tis)
        print("out_s = ", self.tfsm_spec.generate_output_projections(self.s0, self.tis))
        print("out_q = ", self.tfsm_impl.generate_output_projections(self.q0, self.tis))
        print("#######")
        return

class RaceFree_TruncatedTree:
    def __init__(self, tfsm, s0, q0, max_len):
        self.tfsm = tfsm
        self.max_len = int(max_len)
        tis = TimedSequence([])
        self.root = Node(tfsm, s0, q0, s0, q0, tis, max_len)
        return

    def derive_bfs_tree(self):
        queue = [self.root]
        while len(queue) > 0:
            curr_node = queue.pop(0)
            tis = curr_node.derive_successors(queue)
            if tis is not None:
                return tis
        return None

    def print_tree(self):
        queue = [self.root]
        while len(queue) > 0:
            curr_node = queue.pop(0)
            curr_node.print_node()
            queue = queue + list(curr_node.children)
        return

class RaceFree_TruncatedTree_for_testing:
    def __init__(self, tfsm_spec, tfsm_impl, s0, q0, max_len):
        self.tfsm_spec = tfsm_spec
        self.tfsm_impl = tfsm_impl
        self.max_len = int(max_len)
        tis = TimedSequence([])
        self.root = Node_for_testing(tfsm_spec, tfsm_impl, s0, q0, s0, q0, tis, max_len)
        return

    def derive_bfs_tree(self):
        queue = [self.root]
        while len(queue) > 0:
            curr_node = queue.pop(0)
            tis = curr_node.derive_successors(queue)
            if tis is not None:
                return tis
        return None

    def print_tree(self):
        queue = [self.root]
        while len(queue) > 0:
            curr_node = queue.pop(0)
            curr_node.print_node()
            queue = queue + list(curr_node.children)
        return