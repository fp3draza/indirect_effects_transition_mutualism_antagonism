# The role of indirect effects in coevolution as mutualism transitions into antagonism

## Structure of repository 

This repository contains the code and data used to perform the analyses described in "The role of indirect effects in coevolution as mutualism transitions into antagonism". The directory contains three subdirectories: 

* `code`
* `data`

The networks contained in the data directory were downloaded from the [web of life](http://www.web-of-life.es/) online repository. All code was run in `R` version 4.0.2.

To replicate analyses, clone this repository and place it in your home directory. All paths specified in the scripts assume that the repository is at the '~/' level. 

### `code` subdirectory

This directory contains the `R` scripts necessary to recreate our results.

* `coevolution_function.R`: implements the coevolution model.
* `fillfun.R`: converts mutualistic interactions to antagonistic.
* `generate_antagonistic_parameter_dataframe.R`: generates the parameters used to simulate coevolution in the antagonistic networks.
* `generate_mutualistic_parameter_dataframe.R`: generates the parameters used to simulate coevolution in the mutualistic networks.
* `join_coevolution_and_network_structure_results.R`: merges the results from the coevolution simulations performed on the antagonistic and mutualistic networks across all three methods of converting interactions. 
* `plots_manuscript.R`: generates the plots included in the manuscript.
* `process_antagonistic_results.R`: processes the raw output from the coevolution model run on antagonistic networks.
* `process_mutualistic_results.R`: processes the raw output from the coevolution model run on mutualistic networks.
* `process_network_structure_results.R`: calculates PCAs summarising network structures.
* `run_simulation_antagonistic.R`: runs the coevolution simulations on the antagonistic networks.
* `run_simulation_mutualistic.R`: runs the coevolution simulations on the mutualistic networks.
* `structures.R`: calculates network properties.
* `summarise_coevolution_network_structure_results_by_scale.R`: summarises results at the network and species scale.
* `supplementary_plots_manuscript.R`: generates the plots included in as supplementary material.

### `data` subdirectory

This directory contains two subdirectories:

* `data/networks`: contains the antagonistic and mutualistic networks, stored in their respective directories, used in the study.
* `data/parameters`: contains the parameters used to simulate coevolution on the antagonistic and mutualistic networks. 

## General description of workflow

1. Measure network properties: After downloading the mutualistic and antagonistic networks from the [web of life](http://www.web-of-life.es/) online repository, we measure network properites by running `structures.R` and using the `MODULAR` software.
2. Run coevolution simulations: We first generate the parameters needed to simulate coevolution by running `generate_mutualistic_parameter_dataframe.R` and `generate_antagonistic_parameter_dataframe.R`. We then simulate coevolution on all mutualistic and antagonistic networks using `run_simulations_mutualistic.R` and `run_simulations_antagonistic.R` respectively. Internally, these scripts use `coevolution_function.R` and `fillfun.R`.
3. Process raw results: We run `process_mutualistic_results.R` and `process_antagonistic_results.R` to process coevolution results. We run `process_network_structure_results.R` to process network structure results.
4. Data wrangling: We run `join_coevolution_and_network_structure_results.R` to merge all results. We then process the merged results by running `summarise_coevolution_network_structure_results_by_scale.R`.
5. Visualize results: We generate the plots included in the manuscript by running `plots_manuscript.R` and `supplementary_plots_manuscript.R`. Note that these scripts will not write out files by default, they should be stored 'by hand'.
