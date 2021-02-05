# PanGEM Toolbox
The PanGEM Toolbox is intended to aid in the reconstruction of a pangenome scale metabolic network model, and in the retrieval of stoichiometrically balanced, feasible reconstructions for each strain therein. PanGEM Toolbox is intended for use with the Microbial Simulation Environment [MSE](https://github.com/jrcasey/mse_AMT), but the SBML models generated can be used with other systems biology tools. The user may need to implement changes and updates if they intend on running this pipeline for a new taxonomic group. Otherwise, the distribution is currently set up to run a demonstration with 649 strains of *Prochlorococcus*. 

The workflow has been broken up into 7 sequential steps (described below). Unfortunately, due to storage and file size limits on GitHub, genomes must be downloaded either (a) from a static repository (see below) or (b) directly from IMG, NCBI, etc. If you intend to build a PanGEM from scratch (e.g., for another taxonomic group), you will need (b), of course. 

Dependencies: [RAVEN 2.0](https://github.com/SysBioChalmers/RAVEN/wiki), [Mosek 9](https://www.mosek.com/downloads/). 

*You will need to sign up for a free academic license for Mosek.*

### Getting Started
1. Clone the mse_AMT repository to a local directory: `git clone https://github.com/jrcasey/mse_AMT`
2. Download [this zipped directory on Zenodo](http://doi.org/10.5281/zenodo.4477905) containing all *Prochlorococcus* genomes used in this tutorial to `data/genomes/IMG_Dump`
3. In Matlab, add the repository to the path: `addpath(genpath(~/path/to/PanGEM/` and navigate to the directory `cd ~/path/to/PanGEM/`
4. Follow the sequence of steps in the folder `/StepByStep`.
5. Please contact [me](https://jrcasey.github.io/) with questions or comments! 