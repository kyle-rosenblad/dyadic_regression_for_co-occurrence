This repository contains code and data for all analyses accompanying the manuscript “Disentangling trait and phylogenetic signals in species co-occurrence: a Bayesian dyadic regression approach".

Each subdirectory corresponds to a distinct analysis. All subdirectories contain not only the necessary input data and R scripts, but also the intermediate and final output files, unless otherwise noted below. This means that if you prefer, you can run just the later stages of an analysis yourself, picking up where the intermediate output files leave off, without needing to run everything from scratch.

output_files contains the final output files from each of the other subdirectories (described below).

confounding_example generates the plots in Fig. 1, a conceptual example of confounding among species traits in a community assembly analysis.

compnet_workflow_diagram generates the plot in Fig. 2, a conceptual diagram of the compnet workflow.

compnet_ms_metacom_sim_effectsizes runs the metacommunity simulations, and accompanying data analyses, for testing the association between simulated competition strength and the regression model's effect size estimates for pairwise trait distances. There are many R scripts, whose file names all contain integers at the beginning. In general, if you would like to reproduce the whole analysis from scratch, you need to run all “1” scripts before the “2s”, then the “3s”, etc. However, I have included all intermediate output files, so if you download the repository as-is, you can begin running the scripts at any intermediate stage you would like. If you have the opportunity, it is convenient, but by no means necessary, to run all scripts beginning with the same number in parallel. e.g., you can run the “2s” as a batch, then the “3s”, etc.

compnet_ms_metacom_sim_effectsizes_phy runs similar simulations for the analyses of pairwise phylogenetic distances. This and the previous subdirectory produce Fig. 3.

compnet_ms_metacom_sim_falsedetections runs similar simulations for testing the regression model's false detection rate for a "dummy trait" that plays no role in metacommunity dynamics.

compnet_empirical_study runs a compnet analysis of a subset of Kirk et al’s (2022) fish community and trait data, producing Fig. 4, as well as a sensitivity analysis that tests for similar results in other drainage basins.

compnet_sbc runs the simulation-based calibration analysis. The .stan files do not need to be run directly; the R scripts call them. Each of the 8 model configurations has two scripts: one that runs the analysis, and one that processes and visualizes the results. Within each pair, the scripts must be run in this order. One set of intermediate output files is missing. These are the .RDS files containing the compnet models, as these files are too large for GitHub. You can generate these by running the full analysis yourself, (feasible in, at most, a few hours on most personal machines), or feel free to contact me. If the manuscript accompanying this GitHub repository is accepted for publication, I will convert the GitHub repository to a permanent Zenodo repository, which will contain the intermediate .RDS files. As with the metacommunity simulation analyses (see above), you may need to adjust the number of processor cores specified at the beginning of each script.

The metacommunity simulation analyses may be difficult to run from scratch on a personal machine. If you would like to reproduce the whole analysis from scratch, I recommend running them on something like a university cluster (this is what I did) or a commercial cloud computing platform. You may need to adjust the number of processor cores if you don't have enough RAM and/or enough processor cores. (See comments near top of each script regarding how to do this.) If you would only like to run the “4_processing” scripts that follow the metacommunity simulations and data analyses, you can do this smoothly and quickly on any personal machine.

References: Kirk, M. A., Rahel F. J., & Laughlin D. C. (2022). Environmental filters of freshwater fish community assembly along elevation and latitudinal gradients. Global Ecology and Biogeography, 31, 470–485. https://doi.org/10.1111/geb.13439
