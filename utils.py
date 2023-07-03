import math
from scipy import integrate
from pulp import *


class methods:

    def Qi_percent():
        Qi1 = -168
        Qi2 = 16968
        for i in range(77):
            Qi1 += 168
            Qi2-=168
            print('port 1 cargo shipment in month:', Qi1)
            print('port 2 cargo shipment in month:', Qi2)

    def Zi(Pi, Q2, diOd):
        Z2 = Pi * Q2 * diOd
        return Z2

    def my_ui_func(u):
        return math.exp(-u ** 2)

    def BodRo(f, diOd, Di, EF, LF, Qi, V, phi):
        EL_Ro = 2.547 + 4 + 0.03 + 2 + 1.65
        BodRo = (f * diOd * Di * (EF * LF * EL_Ro + phi * EL_Ro)) / (Qi * V * phi * 1e6)
        return BodRo

    def BodRao(f_tilde, dOs, Di, EF_tilde, LF, Qi, V_tilde, phi, f, dOd, V, EF):
        ELRa = 2.547 + 2.01 + 0.134 + 7.37 + 0.402
        ELRo = 2.547 + 4 + 0.03 + 2 + 1.65
        BodRao = (f_tilde * dOs * Di * (EF_tilde * LF * ELRa + phi * ELRa)) / (Qi * V_tilde * phi) + \
                (f * (dOd - dOs) * Di * (EF * LF * ELRo + phi * ELRo)) / (Qi * V * phi * 1e6)
        return BodRao

    def A(alpha, C, BodRo, BodRao):
        A = alpha * C * (BodRo - BodRao)
        return A

    def B(tau, C, BodRo, BodRao):
        B = tau * C * (BodRo - BodRao)
        return B

    def THi(Oi, Di, j):
        THij = Oi * Di ** j
        return THij

    def TC_i(DXi, xi, eps_Ro, eps_Ra, diOs, diOd):

        TCi1 = DXi * xi * (eps_Ro * (diOd - diOs) + eps_Ra * diOs)
        TCi2 = eps_Ro * DXi * xi * diOd
        TCi = TCi1 - TCi2
        return TCi

    def Ki(Di, Qi, ni):
        Ki = (Di - Qi) / ni
        return Ki

    def HRi(R, Ki):
        HRi = R * Ki
        return HRi

    def TCi(TCi0, TCi, HRi):
        TCi = TCi0 + TCi + HRi
        return TCi

    def Ei(Di, Qi, Hi, j):
        Eij = (Di - Qi) * Hi ** j
        return Eij

    def Pi(Zi, A, B, TCi):
        Pi1 = Zi + A - TCi
        Pi2 = Zi - B - TCi
        Pi = Pi1 - Pi2
        return Pi

    def TGR(Di, Pi, diOd, tau, omega, alpha):
        TGR = 0
        for i in range(1, 3):
            TGR += diOd[i-1] * Di[i-1] * Pi[i-1] * (tau - omega - alpha)
        return TGR

    def TSW(pi_list, ui_func, u_start, u_end):
        TSW = sum(pi_list)
        TSW += integrate.quad(ui_func, u_start, u_end)[0]
        return TSW

    def TER(BodRo, BodRao):
        TER = sum(BodRo) + sum(BodRao)
        return TER

    def max_tgr(L_TSW, L_TER):
        prob = LpProblem("Max TGR", LpMaximize)

        TSW = LpVariable("TSW", lowBound=None)
        TER = LpVariable("TER", lowBound=None)
        TGR = LpVariable("TGR", lowBound=None)

        prob += TGR

        prob += TSW >= L_TSW
        prob += TER <= L_TER

        prob.solve()

        status = LpStatus[prob.status]
        tgr = value(TGR.varValue)
        tsw = value(TSW.varValue)
        ter = value(TER.varValue)
        return status, tgr, tsw, ter


    def max_tsw(L_TGR, L_TER):
        prob = LpProblem("Max TSW", LpMaximize)

        TSW = LpVariable("TSW", lowBound=None)
        TER = LpVariable("TER", lowBound=None)
        TGR = LpVariable("TGR", lowBound=None)

        prob += TSW

        prob += TGR >= L_TGR
        prob += TER <= L_TER

        prob.solve()

        status = LpStatus[prob.status]
        tsw = value(TSW.varValue)
        tgr = value(TGR.varValue)
        ter = value(TER.varValue)
        return status, tsw, tgr, ter

    def min_ter(L_TGR, L_TSW):
        prob = LpProblem("Min TER", LpMinimize)

        TSW = LpVariable("TSW", lowBound=None)
        TER = LpVariable("TER", lowBound=None)
        TGR = LpVariable("TGR", lowBound=None)

        prob += TER

        prob += TGR >= L_TGR
        prob += TSW >= L_TSW

        prob.solve()

        status = LpStatus[prob.status]
        ter = value(TER.varValue)
        tgr = value(TGR.varValue)
        tsw = value(TSW.varValue)
        return status, ter, tgr, tsw

#RL
class rl:

    def find_min_max(triad_list):

        min_triad = triad_list[0]
        max_triad = triad_list[0]

        for triad in triad_list:

            if triad < min_triad:
                min_triad = triad

            if triad > max_triad:
                max_triad = triad

        triad_list = max_tgr, max_tsw, min_ter
        min_triad, max_triad = find_min_max(triad_list)

        print(":", min_triad)
        print(":", max_triad)

        return (min_triad, max_triad)

#CART
class cart:

    def find_min_max(triad_list):

        triad_int_list = [int(triad.replace("T", "")) for triad in triad_list]

        min_triad = triad_int_list[0]
        max_triad = triad_int_list[0]

        for triad in triad_int_list:

            if triad < min_triad:
                min_triad = triad

            if triad > max_triad:
                max_triad = triad

        min_triad_str = "T" + str(min_triad)
        max_triad_str = "T" + str(max_triad)

        triad_list = ['max_tgr', 'max_tsw', 'min_ter']

        min_triad, max_triad = find_min_max(triad_list)

        print(":", min_triad)
        print(":", max_triad)

        return (min_triad_str, max_triad_str)

#C4.5
class c45:

    def entropy(class_counts):

        total_count = sum(class_counts)
        entropy = 0

        for count in class_counts:
            if count == 0:
                continue
            frequency = count / total_count
            entropy += -frequency * math.log2(frequency)

        return entropy

    def split_data(data, attribute_index):

        split_data = {}

        for instance in data:
            attribute_value = instance[attribute_index]
            if attribute_value not in split_data:
                split_data[attribute_value] = []
            split_data[attribute_value].append(instance)

        return split_data

    def choose_attribute(data, attribute_indices):

        class_counts = {}

        for instance in data:
            class_label = instance[-1]
            if class_label not in class_counts:
                class_counts[class_label] = 0
            class_counts[class_label] += 1

        base_entropy = entropy(list(class_counts.values()))

        best_attribute_index = None
        best_information_gain = 0

        for attribute_index in attribute_indices:
            attribute_values = set([instance[attribute_index] for instance in data])
            attribute_entropy = 0

            for attribute_value in attribute_values:
                attribute_data = [instance for instance in data if instance[attribute_index] == attribute_value]
                attribute_frequency = len(attribute_data) / len(data)
                attribute_entropy += attribute_frequency * entropy([count for value, count in
                                                                    class_counts.items() if any(instance[-1] == value
                                                                                                for instance in attribute_data)])

            information_gain = base_entropy - attribute_entropy

            if information_gain > best_information_gain:
                best_attribute_index = attribute_index
                best_information_gain = information_gain

        return best_attribute_index

    def majority_class(class_labels):

        class_counts = {}

        for class_label in class_labels:
            if class_label not in class_counts:
                class_counts[class_label] = 0
            class_counts[class_label] += 1

        return max(class_counts, key=class_counts.get)

    def c45(data, attribute_indices):

        class_labels = [instance[-1] for instance in data]

        if len(set(class_labels)) == 1:
            return class_labels[0]

        if not attribute_indices:
            return majority_class(class_labels)

        chosen_attribute_index = choose_attribute(data, attribute_indices)
        tree = {chosen_attribute_index: {}}
        attribute_indices.remove(chosen_attribute_index)

        for attribute_value, attribute_data in split_data(data, chosen_attribute_index).items():
            subtree = c45(attribute_data, attribute_indices[:])
            tree[chosen_attribute_index][attribute_value] = subtree

        return tree

    def find_min_max(triad_list):

        triad_int_list = [int(triad.replace("T", "")) for triad in triad_list]

        data = [[triad] for triad in triad_int_list]

        for instance in data:
            instance.append(0)

        attribute_indices = [0] 
        tree = c45(data, attribute_indices)

        current_node = tree

        while isinstance(current_node, dict):
            attribute_index = list(current_node.keys())[0]
            subtree_key = str(min(triad_int_list)) if "T" + str(min(triad_int_list)) in current_node[attribute_index] else str(max(triad_int_list))
            current_node = current_node[attribute_index][subtree_key]

        min_triad_str = "T" + str(current_node)
        max_triad_str = "T" + str(current_node)

        return (min_triad_str, max_triad_str)