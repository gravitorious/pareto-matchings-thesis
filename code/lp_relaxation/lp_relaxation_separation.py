import random as pyrandom
import csv
from sage.all import *
from fractions import Fraction

# convert to exact rationals
def to_QQ(v, max_den=1000000, tol=1e-8):
    vf = float(v)
    if abs(vf) <= tol: return QQ(0)
    if abs(vf - 1) <= tol: return QQ(1)
    return QQ(Fraction(vf).limit_denominator(max_den))

#feasibility check - all constraints without coalition free constraint
def check_exact_feasibility(vals_QQ, V, A, H, P, Q):
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
    return True

#separation
def separation(sol, G, eps=1e-8):
    """
    Finds a violated coalition-free constraint.
    Input: sol[(a,h)] = current LP value x^*_ah
    Output:
        C = list of vertices [(a,h), ...] forming a violated cycle
        or None if no violated coalition constraint exists.
    """
    #adding weight to each vertex: w(a,h) = 1 - x^*_ah
    w = {v: 1.0 - float(sol[v]) for v in G.vertices()}
    #build weighted directed graph.
    #w(u,v) =  w[v].
    Gw = DiGraph()
    Gw.add_vertices(G.vertices())
    for u, v in G.edge_iterator(labels=False):
        Gw.add_edge(u, v, w[v])
    #all-pairs shortest path distances
    dist = Gw.distance_all_pairs(by_weight=True)
    best_weight = float("inf")
    best_edge = None
    #for every edge (u,v), we check:
    #weight(u,v) + shortest path from v back to u
    for u, v in Gw.edge_iterator(labels=False):
        if v in dist and u in dist[v]:
            cycle_weight = w[v] + dist[v][u]
            if cycle_weight < best_weight:
                best_weight = cycle_weight
                best_edge = (u, v)
    #no directed cycle exists
    if best_edge is None:
        return None
    #no violated cycle
    if best_weight >= 1.0 - eps:
        return None
    #recover cycle
    u, v = best_edge
    path = Gw.shortest_path(v, u, by_weight=True)
    #cycle is (u, v, ..., u)
    #we take the part (v, ..., u)
    C = [u] + path[:-1]
    return C

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
    #generate random objective functions
    pyrandom.seed(1)
    for trial in range(1, 1001):
        objective_coeffs = {v: pyrandom.randint(-20, 20) for v in V}
        lp.set_objective(sum(QQ(objective_coeffs[v]) * x[v[0], v[1]] for v in V))
        try:
            while(True):
                lp.solve()
                sol = lp.get_values(x)
                C = separation(sol, G)
                if C is None:
                    break
                lp.add_constraint(sum(x[a, h] for a, h in C) <= len(C) - 1)
            #rationalize solution
            lp_vals_QQ = {v: to_QQ(sol[v]) for v in V}
            is_fractional = any(QQ(0) < val < QQ(1) for val in lp_vals_QQ.values())
            if is_fractional:
                if check_exact_feasibility(lp_vals_QQ, V, A, H, P, Q):
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
                    return lp_vals_QQ.items()
        except Exception:
            continue
    print("No fractional vertex was found after 1000 trials.")
    return None

def read_preferences(instance):
    P = {}
    houses = set()
    with open("data/" + instance + ".csv", newline="") as prefcsv:
        reader = csv.DictReader(prefcsv)
        for row in reader:
            applicant = int(row["applicant"])
            prefs = []
            for key, value in row.items():
                if key.startswith("pref") and value != "":
                    house = int(value)
                    prefs.append(house)
                    houses.add(house)
            P[applicant] = prefs
    A = sorted(P.keys())
    H = sorted(houses)
    return A, H, P

def lp_solve(instance):
    A, H, P = read_preferences(instance)
    return find_fractional_vertex(A, H, P)