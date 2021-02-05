%% Pangenome assembly - Prochlorococcus
% This script combines the awkward sequence of scripts to create a database
% of strains, genes, KO's, core vs flex to feed into the Smatrix assembly.

% Instead of building the database straight away, we'll just append to a
% single large table, filling out unmatched genes as we go. So we start
% with all the genes from fastas which we can then associate with strains.
% Then we load the results from the KAAS run. Then we match those genes
% with the same verbose gene names as those with KAAS hits and associate to
% the same corresponding KO's. Then we'll import our manually curated list
% of 'still missing genes' and update the table again. 

% The resulting table will then be used to create the rarefaction curves
% and assign core/flex as a column in the table.


%% OPTIONS
SaveMe = 1; % enter 1 to save to Pro_Assembly2



%% Import directory of FASTA's
fileList_path = 'CBIOMES/Pangenomes/Prochlorococcus/NCBIDump/FileList.txt';
fileList = readtable(fileList_path,'delimiter','\t','ReadVariableNames',false);
fileList.Properties.VariableNames = {'FileName'};

%% Import metadata
orgDatabase_path = 'CBIOMES/Pangenomes/Prochlorococcus/AllGenomes3.txt';
orgDatabase = readtable(orgDatabase_path,'delimiter','\t');
varNames = [{'UID'},{'Organism_Name'},	{'Strain'},	{'CladeID'}, {'BioSample'},...	
    {'BioProject'}, {'Assembly'}, {'Size_Mb'}, {'GCpct'}, {'WGS'}, {'Scaffolds'},... 
    {'Genes'}, {'Proteins'},	{'Release_Date'}, {'Modify_Date'}, {'Level'}, ...
    {'RefSeq_FTP'},{'Downloaded_RefSeq'},	{'GenBank_FTP'}];
orgDatabase.Properties.VariableNames = varNames;
orgDatabase = sortrows(orgDatabase,1);

% match organism number, number of scaffolds, and strain name to metadata
for i = 1:numel(fileList.FileName);
    substr = strrep(fileList.FileName{i},'_protein.faa','');
    crap = strfind(orgDatabase.RefSeq_FTP,substr);
    ind = find(not(cellfun('isempty', crap)));
    fileList.orgNumber(i) = orgDatabase.UID(ind);
    fileList.Scaffolds(i) = orgDatabase.Scaffolds(ind);
    fileList.Strain(i) = orgDatabase.Organism_Name(ind);
end % can always add more fields here (e.g., ecotype)


%% Generate database
% Get a list of gene id's and verbose gene names for each fasta file using fastaread
for i = 1:numel(fileList.FileName);
    filename = strcat('CBIOMES/Pangenomes/Prochlorococcus/NCBIDump/',fileList.FileName{i});
    temp2 = fastaread3Headers(filename,'trimheaders',true);
    GeneID{i} = {temp2.Header};
    temp1 = fastaread3(filename);
    GeneVerbose{i} = {temp1.Header};
    Strain{i} = repmat(fileList.Strain(i),1,numel(GeneID{i}));
end

% Create the table and populate with all gene ID's, verbose names, and strains
Dat = struct;
Dat.GeneID = [GeneID{:}]';
% need to remove the gene id and the strain name first
temp1 = [GeneVerbose{:}]';
for i = 1:numel(temp1)
    temp2 = strrep(temp1{i},'MULTISPECIES: ',''); % in case there's a multispecies label
    begin = regexp(temp2,' ');
    begin = begin(1);
    stop = regexp(temp2,'[Proc');
    if isempty(stop)
        stop = regexp(temp2,'[');
    end
    temp3{i} = temp2(begin+1 : stop - 2);
end
   
Dat.GeneVerbose = [temp3(:)];
Dat.Strain = [Strain{:}]';

%% MATCHING STEP 1
% Import KAAS results
KO_path = 'CBIOMES/Pangenomes/Prochlorococcus/KEGG_Orthologies/Combined_KAAS.ko.txt';
KO_dat = readtable(KO_path,'delimiter','\t','ReadVariableNames',false);

% The KO_dat.Var1 gene ID's should be in the same order as Dat.GeneID, so
% just need to append the KO column to Dat.

Dat.KO = KO_dat.Var2;


%% MATCHING STEP 2
% Match verbose gene names which did not get hits in KAAS to those that
% did.

% pull out the set with KO's
KAAS_hitsInd = find(~cellfun(@isempty,Dat.KO));
KAAS_hits = Dat.GeneVerbose(KAAS_hitsInd);
KAAS_hitsKO = Dat.KO(KAAS_hitsInd);


% pull out the set without KO's
KAAS_missesInd = find(cellfun(@isempty,Dat.KO));
KAAS_misses = Dat.GeneVerbose(KAAS_missesInd);

% now do the match
for i = 1:numel(KAAS_misses);
    matchInd = find(strcmp(KAAS_misses{i},KAAS_hits));
    if ~isempty(matchInd);
        Dat.KO(KAAS_missesInd(i)) = KAAS_hitsKO(matchInd(1));
    end
end


%% MATCHING STEP 3
% Import the 'still missing' genes and update the table
stillMissingPath = 'CBIOMES/Pangenomes/Prochlorococcus/RemainingMissingGenes.csv';
stillMissingDat = readtable(stillMissingPath,'delimiter',',','ReadVariableNames',false);

% match
for i = 1:numel(stillMissingDat.Var1);
    stillMissingInd = find(strcmp(stillMissingDat.Var1{i},Dat.GeneID)); % find index
    for j = 1:numel(stillMissingInd);
        if isempty(Dat.KO{stillMissingInd(j)}); % don't want to overwrite so check if there is already a KO present
            Dat.KO{stillMissingInd(j)} = stillMissingDat.Var3{i};
        end
    end
end



%% SAVE?
if SaveMe;
    Pro_assembly2 = Dat;
    Dat.orgDatabase = orgDatabase;
    Dat.fileList = fileList;
    savefast('CBIOMES/Pangenomes/Prochlorococcus/Pro_Assembly2.mat','Dat')
end


