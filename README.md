This repo contains the tex files and source code of my thesis (BSc with integrated MSc) from the University of West Attica, Department of Informatics and Computer Engineering, based on matching theory and especially Pareto matchings in the house allocation problem. My supervisor was Professor Dr. Magos.

After studying some classical matching problems from an algorithmic as well as a polyhedral view, I studied the linear relaxation of the integer program (provided by Eirinakis, P., D. Magos, and I. Mourtos. "A model for House Allocation." 2015 SCinTE: 19). I first showed that the relaxation is not integral and then I constructed a polynomial separation method by formulating the problem in terms of graph theory. I then proposed some candidate cutting inequalities and proved their validity, one of which yields an upper bound on the dimension of the Pareto optimal matching (POM) polytope (convex hull of the POM vectors). Finally, I constructed some facet-defining inequalities for a specific cyclic family of preferences. I also provide the source code I used to perform my experiments, which prints some important statistics, such as all the POM vectors/solutions using the serial dictatorship algorithm, the H-representation and the dimension of the POM polytope, and some results on the linear relaxation, such as a fractional vertex and the inequalities/equalities of the true polytope violated by it, etc. I used SageMath + Python.

The amount of work left is huge. I believe that the linear relaxation of the initial integer program is a really complex object, and it seems (of course) that it is very different from the true POM polytope. The final goal is to come as close as we can to the H-representation of that polytope (or prove how far we can go). Someone who is interested in this problem can benefit from my work by finding some stronger cutting inequalities or, of course, finding facet-defining inequalities that work on more general inputs. Our source code can also be used for that. For example, someone can try to find inequalities from the H-representation that appear more frequently than others and try to find their connection with the input, or even try to prove them. I hope that my work can introduce more people to matching theory and provide some tools to people who do (or decide to do) research on the house allocation problem (or on similar problems). I hope that in the future we will get closer to the H-representation and learn more about it.

## Citation

If you use the results, algorithms, formulations, figures, experimental code, or other material from this repository in academic work, please cite the corresponding thesis:

> N. Mavrogeneiadis, “Matchings with and without preferences, Pareto matchings: problems, structural and polyhedral analysis and algorithmic solution methods,” Integrated Master’s thesis, Dept. Informatics and Computer Engineering, University of West Attica, Athens, Greece, 2026.

The official version of the thesis is available through the institutional repository of the University of West Attica:

https://polynoe.lib.uniwa.gr/xmlui/handle/11400/13037

## License

The thesis and associated academic material in this repository are licensed under the **Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0)** license.

The source code in the [`code/`](code/) directory is licensed separately under the **MIT License**.

See [`LICENSE`](LICENSE) and [`code/LICENSE`](code/LICENSE) for details.


