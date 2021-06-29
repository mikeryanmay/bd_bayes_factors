# Supplementary Archive for: Mistreating tree models as priors compromises our ability to compare phylogenetic tree models

This archive contains code and data necessary to recreate the results for our manuscript, "Mistreating tree models as priors compromises our ability to compare phylogenetic tree models". We performed simulations and empirical data analyses to show the effect of treating tree models as priors on Bayesian phylogenetic model comparison. We explain the relationship between the various files in this archive and our analyses below.

# `yule_bd`

This directory contains simulated data and analysis scripts to recreate the simulation studies described in our Supplementary Material sections S1.1 and S1.2. The subdirectories contain the following files:
* `src` contains scripts for simulating trees and molecular datasets
    - `simulate_trees.R` simulates 10 trees each with S = 8, 16, 32, 64 species
    - `simulate_data.Rev` simulates molecular datasets with 100 sites under a JC model and strict molecular clock on each simulated tree
* `data` contains simulated molecular datasets
    - `JC/bd/n_S` (where S = 8, 16, 32, 64) contains ten simulated molecular datasets for S species simulated with the above scripts
* modules contains `Rev` code for defining all the models we compared in S1.1 and S1.2 (see the jobs directories for scripts to actually perform analyses)
    - `template_sub.Rev` is a template file used to combine various `.Rev` subscripts into one analysis for our comparisons of substitution models (S1.1)
    - `template.Rev` is a template file used to combine various `.Rev` subscripts into one analysis for our comparisons of tree models (S1.2)
    - `sub_models/X_strict.Rev` defines substitution model X (or reversible-jump between JC and K80 in the case of `RJ_strict.Rev`) and a strict molecular clock model
    - `tree_models/Y.Rev` defines the tree model Y (or reversible-jump between Yule and BD tree models in the case of `RJ.Rev`)
    - `analysis/MCMC.Rev` specifies standard MCMC under the defined model (for RJ analysis, also computes the BF between the two models)
    - `analysis/ML.Rev` specifies power-posterior MCMC to estimate the marginal likelihood under the defined model
* `jobs_sub_model_S` (S is the number of species) contains complete `.Rev` scripts for performing the `RevBayes` analyses in S1.1 (comparing molecular substitution models)
* `jobs_tree_model_S` (S is the number of species) contains complete `.Rev` scripts for performing the `RevBayes` analyses in S1.2 (comparing tree models)

# `fbd`

This directory contains simulated data and analysis scripts to recreate the simulation study described in our Supplementary Material section S1.3. The subdirectories contain the following files:
* `src` contains scripts for simulating trees and morphological datasets
    - `simulate_trees.R` simulates 10 trees each with 9 extant and three extinct lineages under a time-heterogeneous fossilized birth-death model (we only used four of the resulting taxon datasets in our study)
    - `simulate_data.Rev` simulates 10 binary morphological datasets with 100 characters under a Mk model and strict morphological clock on each simulated taxon dataset
* `data` contains simulated datasets
    - `dataset_X` contains ten simulated morphological datasets for for tree X, as well as the true tree (`tree.nex`) and stratigraphic data for each species (`taxa.tsv`)
* modules contains `Rev` code for defining all the models we compared in S1.3 (see the jobs directories for scripts to actually perform analyses)
    - `template.Rev` is a template file used to combine various `.Rev` subscripts into one analysis for our comparisons of tree models (S1.2)
    - `sub_models/Mk.Rev` defines a homogeneous Mk model for morphological evolution
    - `tree_models/Y.Rev` defines the tree model Y (one- or two-rate model for fossilization rates; or reversible-jump between Yule and BD tree models in the case of `RJ.Rev`)
    - `analysis/MCMC.Rev` specifies standard MCMC under the defined model (for RJ analysis, also computes the BF between the two models)
    - `analysis/ML.Rev` specifies power-posterior MCMC to estimate the marginal likelihood under the defined model
* `jobs_X` contains complete `.Rev` scripts for performing the `RevBayes` analyses in S1.3 (comparing fossilized birth-death tree models) for each of the X taxon datasets

# Marattiales

This directory contains data and analysis scripts to recreate the empirical analysis described in our Supplementary Material section S1.4. The subdirectories contain the following files:
* `src` contains scripts for simulating taxon datasets for posterior-predictive simulation
    - `simulate_div.R` simulates 50,000 taxon datasets under each model produced by the scripts described below; this script sources the remaining scripts in the `src` directory (which do nothing by themselves)
* `data` contains empirical taxon and morphological datasets
    - `epochs.csv` defines the epochs at which rates change in the diversification models (described below)
    - `ingroup/taxa.tsv` is taxon age data for the ingroup taxa
    - `ingroup/data.nex` is all the binary morphological data for the ingroup taxa
* `modules` contains Rev code for defining all the models we compared in S1.4 (see the jobs directories for scripts to actually perform analyses)
    - `template_empirical.Rev` is a template file used to combine various `.Rev` subscripts into one analysis for our comparisons of tree models (S1.2)
    - `sub_models/Mk_G_relaxed.Rev` defines the Mk model with gamma-distributed rate variation and a relaxed morphological clock (UCLN) model
    - `tree_models/` contains scripts that define the three tree models we used: `constant_rate.Rev` for a model with constant fossilization rates, `phi_variable.Rev` for a model with variable fossilization rates, and `RJ.Rev` for reversible-jump between those two models.
    - `analysis/MCMC.Rev` specifies standard MCMC under the defined model (for RJ analysis, also computes the BF between the two models)
    - `analysis/ML.Rev` specifies power-posterior MCMC to estimate the marginal likelihood under the defined model
* `jobs` contains complete `.Rev` scripts for performing the `RevBayes` analyses in S1.4