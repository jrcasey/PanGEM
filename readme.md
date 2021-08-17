# PanGEM Toolbox
The PanGEM Toolbox is intended to aid in the reconstruction of a pangenome scale metabolic network model, and in the retrieval of stoichiometrically balanced, feasible reconstructions for each strain therein. PanGEM Toolbox is intended for use with the Microbial Simulation Environment [MSE](https://github.com/jrcasey/mse_AMT), but the SBML models generated can be used with other systems biology tools. The user may need to implement changes and updates if they intend on running this pipeline for a new taxonomic group. Otherwise, the distribution is currently set up to run a demonstration with 649 strains of *Prochlorococcus*. 

The workflow has been broken up into 7 sequential steps (described below). Unfortunately, due to file size limits on GitHub (<100MB), several files did not make it into the demo distribution. They are:

1. chem_xref.txt - a link identifier database from MetaNetX, used for annotating metabolites. This can be downloaded from [here](https://www.metanetx.org/mnxdoc/mnxref.html).
2. Assembly\_L1\_[YYYYMMDD] - The level 1 assembly generated by PanGEM Toolbox in Step 2. You will generate this file locally anyway (takes about 10 minutes on my machine), but would have been nice to include in the demo... sigh.
3. IMG_Dump - Genomes must be downloaded either (a) from a static repository (see below) or (b) directly from IMG, NCBI, etc. If you intend to build a PanGEM from scratch (e.g., for another taxonomic group), you will need (b), of course. Download [this zipped directory on Zenodo](http://doi.org/10.5281/zenodo.4477905) containing all *Prochlorococcus* genomes used in this tutorial to `data/genomes/IMG_Dump`

Dependencies: [RAVEN 2.0](https://github.com/SysBioChalmers/RAVEN/wiki)\*, [Mosek 9](https://www.mosek.com/downloads/)\**. 

*\*You will need to sign up for a free academic license for Mosek. Not painless, unfortunately. If you are so motivated, I'd love to hear feedback on using other solvers with this toolbox, especially wrt the performance of solveCS.m.*

*\*\* Very few RAVEN scripts are used, so in a future release we will not require this dependency.*


### Getting Started
1. Clone the mse_AMT repository to a local directory: `git clone https://github.com/jrcasey/mse_AMT`
2. In Matlab, add the repository to the path: `addpath(genpath(~/path/to/PanGEM/` and navigate to the directory `cd ~/path/to/PanGEM/`
3. Follow the sequence of steps in the folder `/StepByStep`.
4. Please contact [me](https://jrcasey.github.io/) with questions or comments! 


### Updates
7/31/2021 - (Under construction) Added proteome synthesis to strain GEMs. This is currently implemented as a final (optional) step in the full PanGEM workflow, as it requires the level 1 assembly containing all the protein sequences. The new Step 8 of the PanGEM workflow generates synthesis reactions for each gene product included in each strain GEM. This is done by parsing the amino acid sequence of each gene product, stored in the level 1 assembly. proteomeGEM's are saved for each strain, along with some information about the stoichiometry of each protein. 