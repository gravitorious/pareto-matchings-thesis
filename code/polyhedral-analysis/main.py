import csv
from sage.all import *
from sd_alg import *

def read_vector_csv(instance):
    path = Path("../results/" + instance + "/" + instance + "_matching_vectors.csv")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, newline="") as vectors:
        reader = csv.DictReader(vectors)
        columns = reader.fieldnames
        #X11 means applicant 1, house 1
        pairs = []
        for col in columns:
            col = col.strip()
            if not col.startswith("X"):
                raise ValueError(f"Invalid column name: {col}")
            a = int(col[1])
            h = int(col[2])
            pairs.append((a, h))
        vertices = []
        for row in reader:
            vec = []
            for col in columns:
                value = row[col].strip()
                if value == "":
                    vec.append(QQ(0))
                else:
                    vec.append(QQ(value))
            vertices.append(vector(QQ, vec))

    #print(pairs)
    #print(vertices)
    return pairs, vertices

def variable_name(pair):
    a, h = pair
    return f"X{a}{h}"

def format_term(coef, name):
    coef = QQ(coef)
    if coef == 1:
        return name
    return f"{coef}*{name}"

def format_sum(terms):
    if not terms:
        return "0"
    return " + ".join(terms)

def print_readable_h_description(polytope, pairs, instance):
    equations = []
    inequalities = []
    for hrep in polytope.Hrepresentation():
        coeffs = list(hrep.A())
        const = QQ(hrep.b())
        lhs_terms = []
        rhs_terms = []
        #hrep means: const + coeffs*x >= 0
        #we rewrite it as:
        #positive terms >= negative terms
        if const > 0:
            lhs_terms.append(str(const))
        elif const < 0:
            rhs_terms.append(str(-const))
        for coef, pair in zip(coeffs, pairs):
            coef = QQ(coef)
            name = variable_name(pair)
            if coef > 0:
                lhs_terms.append(format_term(coef, name))
            elif coef < 0:
                rhs_terms.append(format_term(-coef, name))
        lhs = format_sum(lhs_terms)
        rhs = format_sum(rhs_terms)
        if hrep.is_equation():
            equations.append(f"{lhs} = {rhs}")
        else:
            #change  1 >= X42 + X43 + X45 to X42 + X43 + X45 <= 1
            #if lhs is a/b or -a
            if lhs.replace("/", "").replace("-", "").isdigit() and rhs != "0":
                inequalities.append(f"{rhs} <= {lhs}")
            else:
                inequalities.append(f"{lhs} >= {rhs}")
    path = Path("../results/" + instance + "/" + instance + "_polyhedral_analysis.csv")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", encoding="utf-8") as poltxt:
        poltxt.write("Equations:")
        for eq in equations:
            poltxt.write("  " + eq + "\n")
        poltxt.write("\nInequalities:")
        for ineq in inequalities:
            poltxt.write("  " + ineq + "\n")

if __name__ == "__main__":
    instance = "instance2"
    main_sd(instance)
    pairs, vertices = read_vector_csv(instance)
    convex_hull = Polyhedron(vertices=vertices, base_ring=QQ)
    path = Path("../results/" + instance + "/" + instance + "_polyhedral_analysis.csv")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as poltxt:
        poltxt.write(f"Dimension: {convex_hull.dimension()}\n")
        poltxt.write(f"Number of vertices: {len(vertices)}\n")
    print_readable_h_description(convex_hull, pairs, instance)

