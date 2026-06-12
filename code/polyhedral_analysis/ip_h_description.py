import csv
import itertools
from sage.all import *
from pathlib import Path

def read_vector_csv(instance):
    path = Path("results/" + instance + "/" + instance + "_pom_vectors.csv")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, newline="") as vectors:
        reader = csv.DictReader(vectors)
        columns = reader.fieldnames
        #X11 means applicant 1, house 1
        variables = []
        for col in columns:
            col = col.strip()
            if not col.startswith("X"):
                raise ValueError(f"Invalid column name: {col}")
            a = int(col[1])
            h = int(col[2])
            variables.append((a, h))
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
    return variables, vertices

def variable_name(variable):
    a, h = variable
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

def hrep_sides(hrep, variables):
    coeffs = list(hrep.A())
    const = QQ(hrep.b())
    lhs_terms = []
    rhs_terms = []
    lhs_has_var = False
    rhs_has_var = False
    #sage form: const + coeffs*x >= 0
    #rewrite as: positive terms >= negative terms
    if const > 0:
        lhs_terms.append(str(const))
    elif const < 0:
        rhs_terms.append(str(-const))
    for coef, variable in zip(coeffs, variables):
        coef = QQ(coef)
        name = variable_name(variable)
        if coef > 0:
            lhs_terms.append(format_term(coef, name))
            lhs_has_var = True
        elif coef < 0:
            rhs_terms.append(format_term(-coef, name))
            rhs_has_var = True
    lhs = format_sum(lhs_terms)
    rhs = format_sum(rhs_terms)
    return lhs, rhs, lhs_has_var, rhs_has_var

def readable_hrep(hrep, variables):
    lhs, rhs, lhs_has_var, rhs_has_var = hrep_sides(hrep, variables)
    if hrep.is_equation():
        return f"{lhs} = {rhs}"
    #1 >= X42 + X43 + X45 -> X42 + X43 + X45 <= 1
    #if not lhs_has_var and rhs != "0":
    #    return f"{rhs} <= {lhs}"
    return f"{lhs} >= {rhs}"

def print_h_description(polytope, variables, instance):
    equations = []
    inequalities = []
    for hrep in polytope.Hrepresentation():
        readable = readable_hrep(hrep, variables)
        if hrep.is_equation():
            equations.append(readable)
        else:
            inequalities.append(readable)
    path = Path("results/" + instance + "/" + instance + "_polyhedral_analysis.txt")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a", encoding="utf-8") as poltxt:
        poltxt.write("\nEquations:\n")
        for eq in equations:
            poltxt.write("  " + eq + "\n")
        poltxt.write("\nInequalities:\n")
        for ineq in inequalities:
            poltxt.write("  " + ineq + "\n")

def polyhedral_analysis(instance, variables, vertices):
    convex_hull = Polyhedron(vertices=vertices, base_ring=QQ)
    path = Path("results/" + instance + "/" + instance + "_polyhedral_analysis.txt")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as poltxt:
        poltxt.write(f"Dimension: {convex_hull.dimension()}\n")
        poltxt.write(f"Number of vertices: {len(vertices)}\n")
    print_h_description(convex_hull, variables, instance)
    return convex_hull

#fractional_vertex = dict_items([((1, 1), 1/7), ((1, 2), 6/7)...
def check_fractional_vertex(
    instance,
    convex_hull,
    variables,
    fractional_vertex,
    tol=QQ(0)
):
    """
    checks a fractional relaxation vertex against the H-description
    of the real POM convex hull.
    """
    fractional_vertex = dict(fractional_vertex)
    #build vector in the same order as variables
    #the polytope coordinates follow the variable order stored in `variables`,
    #which is the same order as the columns of the pom-vectors csv.
    #the solver output may be unordered and may omit zero-valued variables.
    #therefore, we build x_vec by iterating over `variables` and retrieving
    #each value from the fractional solution dictionary. the missing entries are 0.
    x_vec = []
    for variable in variables:
        val = fractional_vertex.get(variable, QQ(0))
        x_vec.append(QQ(val))
    x_vec = vector(QQ, x_vec)
    violated_inequalities = []
    satisfied_inequalities = []
    violated_equations = []
    satisfied_equations = []
    for hrep in convex_hull.Hrepresentation():
        #evaluate each H-description constraint at the fractional relaxation vertex.
        #negative value means a violated facet inequality.
        coeffs = vector(QQ, list(hrep.A()))
        const = QQ(hrep.b())
        value = const + coeffs.dot_product(x_vec)
        readable = readable_hrep(hrep, variables)
        if hrep.is_equation():
            if value != 0:
                violated_equations.append((readable, value))
            else:
                satisfied_equations.append((readable, value))
        else:
            if value < tol:
                violated_inequalities.append((readable, value))
            else:
                satisfied_inequalities.append((readable, value))
    path = Path("results/" + instance + "/" + instance + "_fractional_vertex_hdescription_check.txt")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("Fractional relaxation vertex\n")
        f.write("=" * 80 + "\n")
        for variable, val in zip(variables, x_vec):
            if val != 0:
                f.write(f"{variable_name(variable)} = {val}\n")
        f.write("\n")
        f.write("Summary\n")
        f.write("=" * 80 + "\n")
        f.write(f"Number of variables: {len(variables)}\n")
        f.write(f"Number of violated equations: {len(violated_equations)}\n")
        f.write(f"Number of satisfied equations: {len(satisfied_equations)}\n")
        f.write(f"Number of violated inequalities: {len(violated_inequalities)}\n")
        f.write(f"Number of satisfied inequalities: {len(satisfied_inequalities)}\n")
        f.write("\n")
        f.write("Violated equations\n")
        f.write("=" * 80 + "\n")
        for readable, value in violated_equations:
            f.write(f"{readable}\n")
            f.write(f"lhs - rhs = {value}\n\n")
        f.write("\n")
        f.write("Satisfied equations\n")
        f.write("=" * 80 + "\n")
        for readable, value in satisfied_equations:
            f.write(f"{readable}\n")
        f.write("\n")
        f.write("Violated inequalities\n")
        f.write("=" * 80 + "\n")
        for readable, value in violated_inequalities:
            f.write(f"{readable}\n")
            f.write(f"lhs - rhs = {value}\n\n")
        f.write("\n")
        f.write("Satisfied inequalities\n")
        f.write("=" * 80 + "\n")
        for readable, value in satisfied_inequalities:
            f.write(f"{readable}\n")
            f.write(f"lhs - rhs = {value}\n\n")