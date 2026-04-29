import random as pyrandom
from sage.all import *
from fractions import Fraction

def find_fractional_vertex(A, H, P):
    Q = {h: [a for a in A if h in P.get(a, [])] for h in H}
    V = [(a, h) for a in A for h in P.get(a, [])]
    #create envy graph
    G = DiGraph()
    G.add_vertices(V)
    for a1, h1 in V:
        for a2, h2 in V:
            if a1 != a2 and h2 in P.get(a1, []):
                if P[a1].index(h2) < P[a1].index(h1):
                    G.add_edge((a1, h1), (a2, h2))
    raw_cycles = G.all_simple_cycles()
    seen = set()
    plausible_coalitions = []
    for C in raw_cycles:
        C_unique = list(set(C)) #remove the repeated vertex
        if len(C_unique) < 2: continue
        agents = [a for a, h in C_unique]
        houses = [h for a, h in C_unique]
        if len(set(agents)) == len(agents) and len(set(houses)) == len(houses):
            key = tuple(sorted(C_unique))
            if key not in seen:
                seen.add(key)
                plausible_coalitions.append(list(key))
    #lp model
    lp = MixedIntegerLinearProgram(maximization=True, solver="GLPK")
    x = lp.new_variable(real=True, nonnegative=True)
    for a, h in V: lp.set_max(x[a, h], 1.0)
    #add constraints
    for a in A:
        if P.get(a): lp.add_constraint(sum(x[a, h] for h in P[a]) <= 1)
    for h in H:
        if Q[h]: lp.add_constraint(sum(x[a, h] for a in Q[h]) <= 1)
    for a in A:
        for h in P.get(a, []):
            lp.add_constraint(sum(x[a, h_p] for h_p in P[a]) +
                              sum(x[a_p, h] for a_p in Q[h] if a_p != a) >= 1)
    for a in A:
        for h in P.get(a, []):
            worse_houses = [h_p for h_p in P[a] if P[a].index(h) < P[a].index(h_p)]
            lp.add_constraint(sum(x[a_p, h] for a_p in Q[h]) -
                              sum(x[a, h_p] for h_p in worse_houses) >= 0)
    for C in plausible_coalitions:
        lp.add_constraint(sum(x[a, h] for a, h in C) <= len(C) - 1)
    #convert to exact rationals
    def to_QQ(v, max_den=1000000, tol=1e-8):
        vf = float(v)
        if abs(vf) <= tol: return QQ(0)
        if abs(vf - 1) <= tol: return QQ(1)
        return QQ(Fraction(vf).limit_denominator(max_den))
    #feasibility check
    def check_exact_feasibility(vals_QQ):
        for v_item in V:
            if vals_QQ[v_item] < QQ(0) or vals_QQ[v_item] > QQ(1): return False
        for a in A:
            if P.get(a) and sum(vals_QQ[(a, h)] for h in P[a]) > QQ(1): return False
        for h in H:
            if Q[h] and sum(vals_QQ[(a, h)] for a in Q[h]) > QQ(1): return False
        for a in A:
            for h in P.get(a, []):
                s1 = sum(vals_QQ[(a, h_p)] for h_p in P[a])
                s2 = sum(vals_QQ[(a_p, h)] for a_p in Q[h] if a_p != a)
                if s1 + s2 < QQ(1): return False
        for a in A:
            for h in P.get(a, []):
                s1 = sum(vals_QQ[(a_p, h)] for a_p in Q[h])
                worse_houses = [h_p for h_p in P[a] if P[a].index(h) < P[a].index(h_p)]
                s2 = sum(vals_QQ[(a, h_p)] for h_p in worse_houses)
                if s1 - s2 < QQ(0): return False
        for C in plausible_coalitions:
            if sum(vals_QQ[v_item] for v_item in C) > QQ(len(C) - 1): return False
        return True
    #generate random objective functions
    pyrandom.seed(424)
    for trial in range(1, 1001):
        objective_coeffs = {v: pyrandom.randint(-20, 20) for v in V}
        lp.set_objective(sum(QQ(objective_coeffs[v]) * x[v[0], v[1]] for v in V))
        try:
            lp.solve()
            sol = lp.get_values(x)
            #rationalize solution
            lp_vals_QQ = {v: to_QQ(sol[v]) for v in V}
            is_fractional = any(QQ(0) < val < QQ(1) for val in lp_vals_QQ.values())
            if is_fractional:
                if check_exact_feasibility(lp_vals_QQ):
                    print("=" * 60)
                    print(f"Found fractional vertex at trial: {trial}")
                    print("-" * 60)
                    print("Objective coefficients:")
                    for v, c in objective_coeffs.items():
                        if c != 0: print(f"{c} * x_{v}")
                    print("\nFractional Vertex Coordinates:")
                    for v, val_qq in lp_vals_QQ.items():
                        if val_qq > QQ(0):
                            print(f"x_{v} = {val_qq}")
                    return
        except Exception:
            continue
    print("No fractional vertex was found after 1000 trials.")

#input
A_input = ['a1', 'a2', 'a3', 'a4', 'a5']
H_input = ['h1', 'h2', 'h3', 'h4', 'h5']
P_input = {
    'a1': ['h1', 'h2', 'h4'],
    'a2': ['h1', 'h3', 'h2'],
    'a3': ['h3', 'h5', 'h1', 'h4'],
    'a4': ['h3', 'h2', 'h5'],
    'a5': ['h2', 'h1']
}

find_fractional_vertex(A_input, H_input, P_input)