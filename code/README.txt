README - Nikolas Mavrogeneiadis
Code to study Pareto-optimal matchings in the House Allocation problem.

The main idea is the following:
- generate all Pareto-optimal matchings (using Serial Dictatorship),
- convert them to POM vectors,
- build their convex hull,
- extract the H-description (inequalities),
- and compare it with the LP relaxation of the model.

The goal is to understand which inequalities are missing from the LP relaxation,
i.e. which facets of the true polytope cut fractional solutions.

Structure:
data/
    Contains the input instances (preference lists).
sd/
    sd_alg.py
    Implementation of Serial Dictatorship.
    Produces all Pareto-optimal matchings and their vector form.
lp_relaxation/
    lp_relaxation_separation.py
    Solves the LP relaxation of the IP formulation.
    Coalition constraints are added via separation (not all at once).
polyhedral_analysis/
    ip_h_description.py
    builds the convex hull of the POM vectors using Sage.
    Outputs dimension, equations and inequalities.
results/
    For each instance, a folder is created with all outputs.
main.py
    Runs the full pipeline for a given instance.
    Specifically:
    1. runs Serial Dictatorship,
    2. reads the POM vectors,
    3. builds the convex hull,
    4. solves the LP relaxation,
    5. extracts a fractional solution,
    6. checks which convex hull inequalities are violated.

To change the instance:
    instance = "instance1"

Outputs (per instance)

- instance*_pom_vectors.csv
    Incidence vectors of the Pareto-optimal matchings.

- instance*_polyhedral_analysis.txt
    Dimension, equations and inequalities of the convex hull.

- instance*_fractional_vertex_hdescription_check.txt
    Takes a fractional vertex from the LP relaxation and checks
    which inequalities of the convex hull it violates.

Our goal is to observe:
- the structure of the convex hull,
- which inequalities appear,
- and which ones are missing from the LP relaxation.
