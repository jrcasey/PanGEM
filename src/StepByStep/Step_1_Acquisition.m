%% Step 1 - Acquisition

%% 1. Generate a table of organisms within some taxonomic resolution (Here,
% we're using the genus Prochlorococcus). Required column headers include:
% FileID - file name without extension
% Include - Boolean vector. True indicates to include the strain in
% generating a strain model in Step 6
% Strain - Unique name for each strain

% If you are planning to use the MSE toolbox, also include the following
% columns:
% Ecotype - This is used for analysis later on, to partition strains into
% clades. Any heirarchical clustering will do, but you'll need to edit the
% ecotype list to reflect your categories.
% GCpct - This is used to assign coefficients to the DNA component of the
% biomass equation.
% Cell_radius - Cell radius in microns. Used in TpOpt algorithm and
% elsewhere.
% Pmax - This is specific to autotrophs, the maximum rate of photosynthesis
% in fg C cell-1 h-1.
% Genome_size - This is also used to assign coefficients to the biomass
% equation.

%% 2. If you are running the demo with Prochlorococcus, download a folder
% called IMG_Dump from Zenodo (http://doi.org/10.5281/zenodo.4477905) to
% the directory PanGEM/data/genomes/. Otherwise, retrieve protein fasta's from NCBI for each organism, place these all
% in a single directory.

%% 3. Concatenate all fasta's into a single file by navigating (in a
% terminal) to the directory where they are located and running cat *.faa*
% > combined.faa 

% I've included a bash script to help with unpacking tarballs and moving
% things around, in case it helps. 

%% 4. In the same directory, generate a list of all of the fasta filenames
% (not including the concatenated file) in a tab-delimited text file called
% 'FileList.txt' 

%% 5. Submit the combined.faa file to KEGG's KAAS server for annotation.
% This can be done at: https://www.genome.jp/kaas-bin/kaas_main. Use GHOSTX
% mapping and include a list of closely related organisms and preferably
% some that are model organisms. For Prochlorococcus, I used the following
% list of KEGG organism ID's: (pae, cre, eco, noc, mej, bma, hpy, gsu, pub,
% pel, bsu, sau, syz, syy, syw, syc, syf, syh, mpp, cvr, sce, ath, pma,
% pmm, pmt, pmn, pmi, pmb, pmc, pmf, pmg, pmh, pmj, pme, prc, prm, cyt)
% After a while, you will get an email with a link to your annotated
% results file, which will be a tab-delimited list with two columns, the
% first is the NCBI gene ID and the second is the corresponding KEGG
% ortholog (KO). Save this file as 'Combined_KAAS.ko.txt' in the same
% directory (FileNames.FastaDir)
% 

