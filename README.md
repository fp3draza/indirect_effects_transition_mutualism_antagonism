# The role of indirect effects in coevolution along the mutualism-antagonism continuum

## About this respository

This repository contains the code and data used to perform the analyses described in *"The role of indirect effects in coevolution along the mutualism-antagonism continuum"*. The directory contains four subdirectories:

-   `code`
-   `data`
-   `output`
-   `results`

The networks contained in the data directory were downloaded from the [web of life](http://www.web-of-life.es/) online repository. All code was run in `R` version 4.0.2.

To replicate analyses, clone this repository and place it in your home directory. All paths specified in the scripts assume that the repository is at the home level (`'~/'`). Note that the scripts will create further repositories inside this one as needed.

### `code` subdirectory

This directory contains the `R` scripts necessary to recreate our results.

-   `coevolution_function.R`: implements the coevolution model (this is an auxiliary function, used by other scripts).
-   `fillfun.R`: converts mutualistic interactions to antagonistic (this is an auxiliary function, used by other scripts).
-   `fit_linear_models.R`: fits linear model to the results in order to visualise trends (this is an auxiliary function, used by other scripts).
-   `functions_for_trait_matching_figure.R`: returns the output of running the coevolution model once on a given simulation. These results are used to construct a figure (this is an auxiliary function, used by other scripts).
-   `generate_antagonistic_parameter_dataframe.R`: generates the parameters used to simulate coevolution in the antagonistic networks.
-   `generate_mutualistic_parameter_dataframe.R`: generates the parameters used to simulate coevolution in the mutualistic networks.
-   `join_coevolution_and_network_structure_results.R`: merges the results from the coevolution simulations performed on the antagonistic and mutualistic networks and the data on the structural properties of the networks.
-   `measure_network_properties.R`: computes a set of structural descriptors of all mutualistic and antagonistic networks.
-   `plots_manuscript.R`: generates the plots included in main text of the manuscript.
-   `process_antagonistic_results.R`: processes the raw output from the coevolution model run on antagonistic networks.
-   `process_mutualistic_results.R`: processes the raw output from the coevolution model run on mutualistic networks.
-   `process_network_structure_results.R`: merges the results on the structural descriptors of antagonistic and mutualistic networks.
-   `run_simulation_antagonistic.R`: runs the coevolution simulations on the antagonistic networks.
-   `run_simulation_mutualistic.R`: runs the coevolution simulations on the mutualistic networks.
-   `summarise_coevolution_network_structure_results_by_scale.R`: summarises results at the network and species scale.
-   `supplementary_plots_manuscript.R`: generates the plots included as supplementary material to the manuscript.

### `data` subdirectory

This directory contains two subdirectories:

-   `data/networks`: contains the antagonistic and mutualistic networks, stored in their respective directories, used in the study.
-   `data/parameters`: contains the parameters used to simulate coevolution on the antagonistic and mutualistic networks. These files are created after running `generate_antagonistic_parameter_dataframe.R` and `generate_mutualistic_parameter_dataframe.R`.

### `output` subdirectory

This directory contains three subdirectories:

-   `output/mutualistic_to_antagonistic`: contains the output of running the coevolution model on a set of mutualistic networks that are progressively shifted to become antagonistic. This directory contains a subdirectory for each combination of fraction of antagonists in the networks and the strategy used to assign antagonists. For example, the output of converting 20% of interactions in a network to antagonism starting from the most specialist species is stored in the subdirectory: `all_networks_m_07_a_02_specialist/`. This directory and its contents are generated after running `run_simulation_mutualistic.R`.

-   `output/antagonistic_to_mutualistic`: contains the output of running the coevolution model on a set of antagonistic networks that are progressively shifted to become mutualistic. This directory contains a subdirectory for each combination of fraction of antagonists in the networks and the strategy used to assign antagonists. For example, the output of converting 20% of interactions in a network to mutualstic starting from the most specialist species is stored in the subdirectory: `all_networks_m_07_a_08_specialist/`. This directory and its contents are generated after running `run_simulation_antagonistic.R`.

-   `output/network_structure`: contains a pair of files detailing the structural properties of the mutualistic and antagonistic networks. This directory and its contests are generated after running `measure_network_properties.R`.

### `results` subdirectory

This directory contains a set of files that contain the post-processed output of the coevolution model and network structural descriptors. The files are generated after running either `process_mutualistic_results.R`, `process_antagonistic_results.R` or `process_network_structure_results.R`.

## General description of workflow

1.  Measure network properties: After downloading the mutualistic and antagonistic networks from the [web of life](http://www.web-of-life.es/) online repository, we measure network properites by running `measure_network_properties.R`.
2.  Run coevolution simulations: We first generate the parameters needed to simulate coevolution by running `generate_mutualistic_parameter_dataframe.R` and `generate_antagonistic_parameter_dataframe.R`. We then simulate coevolution on all mutualistic and antagonistic networks using `run_simulations_mutualistic.R` and `run_simulations_antagonistic.R` respectively. Internally, these scripts use `coevolution_function.R` and `fillfun.R`.
3.  Process raw results: We run `process_mutualistic_results.R` and `process_antagonistic_results.R` to process coevolution results. We run `process_network_structure_results.R` to process network structure results.
4.  Data wrangling: We run `join_coevolution_and_network_structure_results.R` to merge all results. We then process the merged results by running `summarise_coevolution_network_structure_results_by_scale.R`.
5.  Visualise results: We generate the plots included in the manuscript by running `plots_manuscript.R` and `supplementary_plots_manuscript.R`. Note that these scripts will not store the figures. Internally, these scripts use `fit_linear_models.R` and `functions_for_trait_matching_figure.R` to generate the figures.

## Sensitivity analyses

We provide the data to generate the figures showing the sensitivity analyses shown in the manuscript (figures A10-A12) in `results/sensitivity_analyses`. We do not provide an Rscript to generate this data. Instead, it should be generated using the workflow detailed above and by implementing the following tweaks:

-   To generate the data exploring the effect of the strength of coevolution (figure A10), change the `m_mean` parameter in `generate_mutualistic_parameter_dataframe.R` and `generate_antagonistic_parameter_dataframe.R` and the `m` parameter in `run_simulation_mutualistic` and `run_simulation_antagonistic`.

-   To generate the data exploring the effect of the sensitivity to trait matching (figure A11), change the `alpha` parameter in `run_simulation_mutualistic` and `run_simulation_antagonistic`.

-   To generate the data exploring the effect of the sensitivity of resource species to consumer species (figure A12), change the `epsilon` parameter in `run_simulation_mutualistic` and `run_simulation_antagonistic`.
