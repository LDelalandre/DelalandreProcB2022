.
├── data 
|	├── code_ForCEEPS_simulations # Code to perform simulations from the ForCEEPS model (hosted on Capsis: http://capsis.cirad.fr/capsis/help_en/forceeps, contact them for any informations about the platform)
|	├── raw # Raw data. /!\ Contains raw simulation outputs, which are not hosted on the GitHub repository. Script to download them in "final/1.2. Import raw data.R"
|	├── processed # Processed data. Available from GitHub. /!\ Only generated from: "final/2. Process data.R" if raw data have been downloaded and de-zipped. 
|
├── docs                    # Documentation files (alternatively `doc`)
├── src                     # Source files (alternatively `lib` or `app`)
├── test                    # Automated tests (alternatively `spec` or `tests`)
├── tools                   # Tools and utilities
├── LICENSE.txt
└── README.txt


A dire dans le README:

Code for simulations forest dynamics is to be found in folder
Output of ForCEEPS simulations

To perform the analysis:

\Final
* 1.Distinctiveness computation_command files.R
i) Computes the distinctiveness of all the species. Writes it in a .txt file.
ii) Writes the command files for the four sites, removing species in increasing, decreasing, and random order of distinctiveness.
NB: seeds are set from 1 to 10 for the ten random orders generated.

* 2.Process_data.R



\R
* Before simulations.R 
1. Computes the distinctiveness of all the species. Writes it in a .txt file.
2. Writes the command files for the four sites, removing species in increasing, decreasing, and random order of distinctiveness.

* Analysis_data.R

* Monocultures_functions.R

\Analysis
