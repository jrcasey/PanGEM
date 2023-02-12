%% Step 6c - Random KO Models 
% Generate Random KO selected GEMs from PanGEM (both standalone format and PanGEM
% format)
% Loads current (annotated) PanGEM and runs CS algorithm to extract
% strain-specific models.
% You can edit solveCS to be less conservative than I have done here. The
% gVec vector of biomass formation rate lower bounds will index the
% reduction in growth rate that you're looking for. Currently, I have it
% set to search for the sparse set across the whole gVec, but one could
% argue...

% As an optional step, each strain model is tested for feasibility under a
% number of different limitations (nitrogen, light, etc). 

% (Less than 1 minute runtime)

% Tatsuro Tanioka 20220725- Modified "Step_6_Strain_Models.m" by John Casey

clear

%%
totalKOs = [910, 920, 930, 940, 950];

% Selecting total random KO profiles generated
totalrndsamples = 100;
totalKO_names = strcat("Random_top",string(totalKOs));
%% Paths
working_dir  = '/Users/tatsurotanioka/Desktop/Project/BioGOSHIP_MetaPanGEM/';
FileNames.PanGEM_path = append(working_dir, 'PanGEM/data/models/Pro_PanGEM.mat');
FileNames.Destination_RandomMod = append(working_dir, 'Metagenome_GEM/models/RandomMod_',num2str(totalrndsamples),'.mat');
FileNames.Destination_PAgenes_RandomMod = append(working_dir, 'Metagenome_GEM/models/PAgenes_RandomMod_',num2str(totalrndsamples),'.mat');
%% Create RandomMod for different throsholds
FileNames.L2Assembly_path = append(working_dir, 'Make_KO_Table/KO_PresAbs_Table/MATfiles/','Random_Assembly.mat');
FileNames.Destination_L3 = append(working_dir, 'Metagenome_GEM/Assembly_L3/','Ramdom_Assembly_L3','.mat');

%% add PanGEM codes to the path
addpath(genpath(append(working_dir,'PanGEM/src/')))

%% Import PanGEM and the L2 Assembly
load(FileNames.PanGEM_path)
load(FileNames.L2Assembly_path)

sample = Meta_Assembly.SampleID;
nAlpha = numel(totalKO_names);
Meta_Assembly.uniqueKO = cellstr(Meta_Assembly.uniqueKO);
%% Load essential reaction subset from PanGEM
% Precauclated in Step6_Strain_Models.m so use the same list
load(append(working_dir,'PanGEM/data/models/essentialRxns.mat'));

%% Get reactions to remove for each metagenome sample and for each confidence level
for a = 1:totalrndsamples
    for b = 1:nAlpha
        [rxnsToRemove{a,b}] = getRxnsToRemove_metagenomes(PanGEM,Meta_Assembly, sample{a}, totalKO_names(b));
    end
end

%% Exclude essential reactions from reactions to remove and generate metagenome based models
load(append(working_dir,'PanGEM/data/models/essentialRxns.mat'));

for a = 1:totalrndsamples
    for b = 1:nAlpha
        tempRxnstoRemove = rxnsToRemove{a,b};
        canRemove{a,b} = setdiff(tempRxnstoRemove,PanGEM.rxns(essentialRxns));
    end
end
PAMat_model = zeros(numel(PanGEM.genes),totalrndsamples,nAlpha); 
for a = 1:totalrndsamples
    for b = 1:nAlpha
        targSample = sample{a};
        topSample = totalKO_names{b};
        [junk, removeIdx] = intersect(PanGEM.rxns,canRemove{a,b});
        % RandomMod(b).(targSample) = PanGEM;
        RandomMod.(topSample).(targSample) = PanGEM;
        % zero out S matrix columns corresponding to these reactions
        % RandomMod(b).(targSample).S(:,removeIdx) = 0;
        RandomMod.(topSample).(targSample).S(:,removeIdx) = 0;
        % store geneset in a special model field (this is to avoid indexing
        % issues with actually removing genes)
        % sampleRxns_Idx = find(sum(abs(RandomMod(b).(targSample).S),1));
        sampleRxns_Idx = find(sum(abs(RandomMod.(topSample).(targSample).S),1));
        [junk, sampleGenes_Idx] = find(PanGEM.rxnGeneMat(sampleRxns_Idx,:));
        % RandomMod(b).(targSample).geneSubset = unique(PanGEM.genes(sampleGenes_Idx));
        RandomMod.(topSample).(targSample).geneSubset = unique(PanGEM.genes(sampleGenes_Idx));

        % tempSol = solveLP(RandomMod(b).(targSample));
        tempSol = solveLP(RandomMod.(topSample).(targSample));
        if tempSol.stat
            mu1(a,b) = -tempSol.f;
        end
        % [junk, include_idx] = intersect(PanGEM.genes,RandomMod(b).(targSample).geneSubset);
        [junk, include_idx] = intersect(PanGEM.genes,RandomMod.(topSample).(targSample).geneSubset);
        PAMat_model(include_idx,a,b) = 1;
        PAgenes_RandomMod(b).PresenceAbsenceMatrix = PAMat_model(:,:,b);
    end
end

%% 
modelFolder = append(working_dir, 'Metagenome_GEM/models/','RandomMod_standalone');
if ~exist(modelFolder, 'dir')
    mkdir(modelFolder)
end

for b = 1:nAlpha
    for a = 1:totalrndsamples
        targSample = sample{a};
        % Get indices of genes and reactions to remove. 
        [junk, removeGenes_Idx] = setdiff(PanGEM.genes, RandomMod.(topSample).(targSample).geneSubset);
        [junk, removeRxn_Idx] = intersect(PanGEM.rxns,canRemove{a,b});

        % Initialize metagenome model
        RandomMod_Standalone.(targSample) = PanGEM;

        % Gene fields
        RandomMod_Standalone.(targSample).genes(removeGenes_Idx) = [];
        RandomMod_Standalone.(targSample).rxnGeneMat(:,removeGenes_Idx) = [];
        RandomMod_Standalone.(targSample).geneMiriams(removeGenes_Idx) = [];
        RandomMod_Standalone.(targSample).geneSBOTerms(removeGenes_Idx) = [];

        % Reaction fields
        RandomMod_Standalone.(targSample).S(:,removeRxn_Idx) = [];
        RandomMod_Standalone.(targSample).rxns(removeRxn_Idx) = [];
        RandomMod_Standalone.(targSample).lb(removeRxn_Idx) = [];
        RandomMod_Standalone.(targSample).ub(removeRxn_Idx) = [];
        RandomMod_Standalone.(targSample).c(removeRxn_Idx) = [];
        RandomMod_Standalone.(targSample).rules(removeRxn_Idx) = [];
        RandomMod_Standalone.(targSample).grRules(removeRxn_Idx) = [];
        RandomMod_Standalone.(targSample).rxnGeneMat(removeRxn_Idx,:) = [];
        RandomMod_Standalone.(targSample).rxnNames(removeRxn_Idx) = [];
        RandomMod_Standalone.(targSample).subSystems(removeRxn_Idx) = [];
        RandomMod_Standalone.(targSample).eccodes(removeRxn_Idx) = [];
        RandomMod_Standalone.(targSample).rxnComps(removeRxn_Idx) = [];
        RandomMod_Standalone.(targSample).rxnKEGGID(removeRxn_Idx) = [];
        RandomMod_Standalone.(targSample).rxnMetaNetXID(removeRxn_Idx) = [];
        RandomMod_Standalone.(targSample).rxnSBOTerms(removeRxn_Idx) = [];
        RandomMod_Standalone.(targSample).rev(removeRxn_Idx) = [];

        RandomMod_Standalone.(targSample).rxnKEGGID = RandomMod_Standalone.(targSample).rxnKEGGID';
        % Remove unused metabolites
        unusedMet_Idx = find(sum(abs(RandomMod_Standalone.(targSample).S),2)==0);

        % Metabolite fields
        RandomMod_Standalone.(targSample).S(unusedMet_Idx,:) = [];
        RandomMod_Standalone.(targSample).mets(unusedMet_Idx) = [];
        RandomMod_Standalone.(targSample).b(unusedMet_Idx) = [];
        RandomMod_Standalone.(targSample).csense(unusedMet_Idx) = [];
        RandomMod_Standalone.(targSample).metCharges(unusedMet_Idx) = [];
        RandomMod_Standalone.(targSample).metFormulas(unusedMet_Idx) = [];
        RandomMod_Standalone.(targSample).metNames(unusedMet_Idx) = [];
        RandomMod_Standalone.(targSample).metComps(unusedMet_Idx) = [];    
        RandomMod_Standalone.(targSample).metKEGGID(unusedMet_Idx) = [];
        RandomMod_Standalone.(targSample).metPubChemID(unusedMet_Idx) = [];
        RandomMod_Standalone.(targSample).metChEBIID(unusedMet_Idx) = [];
        RandomMod_Standalone.(targSample).metInChI(unusedMet_Idx) = [];
        RandomMod_Standalone.(targSample).metMetaNetXID(unusedMet_Idx) = [];
        RandomMod_Standalone.(targSample).metSBOTerms(unusedMet_Idx) = [];

        % Annotate model field
        RandomMod_Standalone.(targSample).description = targSample;
        RandomMod_Standalone.(targSample).modelID = strcat('iTT',mat2str(numel(RandomMod.(topSample).(targSample).geneSubset)),'_',targSample,'_',totalKO_names{b});
        RandomMod_Standalone.(targSample).id = strcat(targSample,'_',totalKO_names{b});

       % Export standalone model to mat
        outputFile = append(modelFolder,'/',targSample,'_',totalKO_names{b},'.mat');
        model = RandomMod_Standalone.(targSample);
        mod = model;
        % save(outputFile,'model');
        save(outputFile,'mod');
    end
end

%% Get PAMat for reactions
PAMat_rxns = ones(numel(PanGEM.rxns),totalrndsamples,nAlpha);
for a = 1:totalrndsamples
    for b = 1:nAlpha
        [junk, junk2, zeroOut_Idx] = intersect(canRemove{a,b},PanGEM.rxns);
        PAMat_rxns(zeroOut_Idx,a,b) = 0;
    end
end
%% Save metagenome sample models and L3 structure
Meta_Assembly_L3 = Meta_Assembly;
Meta_Assembly_L3.PAMat_rxns = PAMat_rxns;
%%
save(FileNames.Destination_L3,'Meta_Assembly_L3');
save(FileNames.Destination_RandomMod,'RandomMod');
save(FileNames.Destination_PAgenes_RandomMod,'PAgenes_RandomMod');














