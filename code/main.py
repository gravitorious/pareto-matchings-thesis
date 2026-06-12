import polyhedral_analysis.ip_h_description as ip_convex_hull
import sd.sd_alg as sd
import lp_relaxation.lp_relaxation_separation as relaxation

if __name__ == "__main__":
    instance = "instance1"
    sd.sd_alg(instance) #run serial dictatorship algorithm
    variables, vertices = ip_convex_hull.read_vector_csv(instance) #read POM vectors
    convex_hull = ip_convex_hull.polyhedral_analysis(instance, variables, vertices)
    fractional_solution = relaxation.lp_solve(instance)
    if fractional_solution is None:
        print("No fractional vertex was found")
    else:
        ip_convex_hull.check_fractional_vertex(instance, convex_hull, variables, fractional_solution)
