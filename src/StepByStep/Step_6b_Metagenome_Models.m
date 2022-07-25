%% Step 6b - Metagenome Models 
% Generate Metagenome GEMs from PanGEM (both standalone format and PanGEM
% format)
% Loads current (annotated) PanGEM and runs CS algorithm to extract
% metagenome-specific models.
% You can edit solveCS to be less conservative than I have done here. The
% gVec vector of biomass formation rate lower bounds will index the
% reduction in growth rate that you're looking for. Currently, I have it
% set to search for the sparse set across the whole gVec, but one could
% argue...

% As an optional step, each strain model is tested for feasibility under a
% number of different limitations (nitrogen, light, etc). 

% (Less than 1 minute runtime)
% Tatsuro Tanioka 20220725- Modified "Step_6_Strain_Models.m" by John Casey

%% 
clear
% Set transect here
% Transect = "C13"
% Transect = "I09"
% Transect = "I07"
Transect = "Globe_AMT_P18_etc"
%% Paths
working_dir  = '/Users/tatsurotanioka/Desktop/Project/BioGOSHIP_MetaPanGEM/';
FileNames.PanGEM_path = append(working_dir, 'PanGEM/data/models/Pro_PanGEM.mat');
FileNames.L2Assembly_path = append(working_dir, 'Make_KO_Table/KO_PresAbs_Table/MATfiles/', Transect, '_Assembly.mat');
FileNames.Destination_L3 = append(working_dir, 'Metagenome_GEM/',Transect,'_Assembly_L3.mat');
FileNames.Destination_SampleMod = append(working_dir, 'Metagenome_GEM/models/',Transect,'_SampleMod.mat');
FileNames.Destination_SampleMod_Standalone = append(working_dir, 'Metagenome_GEM/models/',Transect,'_SampleMod_Standalone.mat');

% FileNames.PanGEM_path = 'data/models/Pro_PanGEM.mat';
% FileNames.L2Assembly_path = 'data/assemblies/Assembly_L2_20200830.mat';
% FileNames.Destination_L3 = 'data/assemblies/Assembly_L3_20200830.mat';
% FileNames.Destination_StrMod = 'data/models/StrMod.mat';
% FileNames.Destination_StrMod_Standalone = 'data/models/StrMod_Standalone.mat';

%% Import PanGEM and the L2 Assembly
load(FileNames.PanGEM_path)

load(FileNames.L2Assembly_path)
%%
% strain = Pro_Assembly_L2.orgDatabase.StrainName(find(Pro_Assembly_L2.orgDatabase.Include));
% nStrains = numel(strain);

sample = Meta_Assembly.SampleID;
nSample = numel(sample);

%% Get essential reaction subset from PanGEM
% BIOMASSMinimal contains the minimum set of BOF compounds required by all
% strains to produce biomass. Any additional components are considered
% non-essential

% minBOFind = find(strcmp('BIOMASSMinimal',PanGEM.rxns));

% BIOMASSCRUDE contains the full set of BOF compounds required by all
% strains to produce biomass. This will yield the most conservative
% essential subset.
minBOFind = find(strcmp('BIOMASSCRUDE',PanGEM.rxns));

PanGEM.c = zeros(numel(PanGEM.rxns),1);
PanGEM.c(minBOFind) = 1;

[essentialRxns, results] = solveCS(PanGEM);
% also don't mess with photosynthesis, oxphos and gap-filled reactions...
% so let's add those to the list
gapFilled = find(cellfun(@isempty,PanGEM.grRules));
PSInd = find(strcmp('Photosynthesis',PanGEM.subSystems));
OxPhosInd = find(strcmp('Oxidative phosphorylation',PanGEM.subSystems));
essentialRxns2 = unique([essentialRxns, PSInd',OxPhosInd',gapFilled']);
essentialRxns = essentialRxns2;

% specific (e.g., experimentally known) reactions to remove from this list
rmList = [{'R00475'}]; % glycolate oxidase
for a = 1:numel(rmList)
    rmIdx(a) = find(strcmp(rmList,PanGEM.rxns));
end
essentialRxns(find(ismember(essentialRxns,rmIdx))) = [];

%% Get reactions to remove for each metagenome sample
   
% for a = 1:nStrains
%     [rxnsToRemove{a}] = getRxnsToRemove(PanGEM,Pro_Assembly_L2, strain{a});
% end
%%
for a = 1:nSample
    [rxnsToRemove{a}] = getRxnsToRemove_metagenomes(PanGEM,Meta_Assembly, sample{a}, 'alpha_95');
end

%% Exclude essential reactions from reactions to remove and generate metagenome based models
for a = 1:nSample
    tempRxnstoRemove = rxnsToRemove{a};
    canRemove{a} = setdiff(tempRxnstoRemove,PanGEM.rxns(essentialRxns));
end

for a = 1:nSample
    targSample = sample{a};
    [junk, removeIdx] = intersect(PanGEM.rxns,canRemove{a});
    SampleMod.(targSample) = PanGEM;
    % zero out S matrix columns corresponding to these reactions
    SampleMod.(targSample).S(:,removeIdx) = 0;
    % store geneset in a special model field (this is to avoid indexing
    % issues with actually removing genes)
    sampleRxns_Idx = find(sum(abs(SampleMod.(targSample).S),1));
    [junk, sampleGenes_Idx] = find(PanGEM.rxnGeneMat(sampleRxns_Idx,:));
    SampleMod.(targSample).geneSubset = unique(PanGEM.genes(sampleGenes_Idx));
    
    tempSol = solveLP(SampleMod.(targSample));
    if tempSol.stat
        mu1(a) = -tempSol.f;
    end
end

%% 
for a = 1:nSample
    targSample = sample{a};
    % Get indices of genes and reactions to remove. 
    [junk, removeGenes_Idx] = setdiff(PanGEM.genes, SampleMod.(targSample).geneSubset);
    [junk, removeRxn_Idx] = intersect(PanGEM.rxns,canRemove{a});
    
    % Initialize metagenome model
    SampleMod_Standalone.(targSample) = PanGEM;
    
    % Gene fields
    SampleMod_Standalone.(targSample).genes(removeGenes_Idx) = [];
    SampleMod_Standalone.(targSample).rxnGeneMat(:,removeGenes_Idx) = [];
    SampleMod_Standalone.(targSample).geneMiriams(removeGenes_Idx) = [];
    SampleMod_Standalone.(targSample).geneSBOTerms(removeGenes_Idx) = [];
    
    % Reaction fields
    SampleMod_Standalone.(targSample).S(:,removeRxn_Idx) = [];
    SampleMod_Standalone.(targSample).rxns(removeRxn_Idx) = [];
    SampleMod_Standalone.(targSample).lb(removeRxn_Idx) = [];
    SampleMod_Standalone.(targSample).ub(removeRxn_Idx) = [];
    SampleMod_Standalone.(targSample).c(removeRxn_Idx) = [];
    SampleMod_Standalone.(targSample).rules(removeRxn_Idx) = [];
    SampleMod_Standalone.(targSample).grRules(removeRxn_Idx) = [];
    SampleMod_Standalone.(targSample).rxnGeneMat(removeRxn_Idx,:) = [];
    SampleMod_Standalone.(targSample).rxnNames(removeRxn_Idx) = [];
    SampleMod_Standalone.(targSample).subSystems(removeRxn_Idx) = [];
    SampleMod_Standalone.(targSample).eccodes(removeRxn_Idx) = [];
    SampleMod_Standalone.(targSample).rxnComps(removeRxn_Idx) = [];
    SampleMod_Standalone.(targSample).rxnKEGGID(removeRxn_Idx) = [];
    SampleMod_Standalone.(targSample).rxnMetaNetXID(removeRxn_Idx) = [];
    SampleMod_Standalone.(targSample).rxnSBOTerms(removeRxn_Idx) = [];
    SampleMod_Standalone.(targSample).rev(removeRxn_Idx) = [];
    
    SampleMod_Standalone.(targSample).rxnKEGGID = SampleMod_Standalone.(targSample).rxnKEGGID';
    % Remove unused metabolites
    unusedMet_Idx = find(sum(abs(SampleMod_Standalone.(targSample).S),2)==0);

    % Metabolite fields
    SampleMod_Standalone.(targSample).S(unusedMet_Idx,:) = [];
    SampleMod_Standalone.(targSample).mets(unusedMet_Idx) = [];
    SampleMod_Standalone.(targSample).b(unusedMet_Idx) = [];
    SampleMod_Standalone.(targSample).csense(unusedMet_Idx) = [];
    SampleMod_Standalone.(targSample).metCharges(unusedMet_Idx) = [];
    SampleMod_Standalone.(targSample).metFormulas(unusedMet_Idx) = [];
    SampleMod_Standalone.(targSample).metNames(unusedMet_Idx) = [];
    SampleMod_Standalone.(targSample).metComps(unusedMet_Idx) = [];    
    SampleMod_Standalone.(targSample).metKEGGID(unusedMet_Idx) = [];
    SampleMod_Standalone.(targSample).metPubChemID(unusedMet_Idx) = [];
    SampleMod_Standalone.(targSample).metChEBIID(unusedMet_Idx) = [];
    SampleMod_Standalone.(targSample).metInChI(unusedMet_Idx) = [];
    SampleMod_Standalone.(targSample).metMetaNetXID(unusedMet_Idx) = [];
    SampleMod_Standalone.(targSample).metSBOTerms(unusedMet_Idx) = [];
    
    % Annotate model field
    SampleMod_Standalone.(targSample).description = targSample;
    SampleMod_Standalone.(targSample).modelID = strcat('iTT',mat2str(numel(SampleMod.(targSample).geneSubset)),'_',targSample);
    SampleMod_Standalone.(targSample).id = targSample;

    % Export standalone model to mat
    modelFolder = append(working_dir, 'Metagenome_GEM/models/',Transect,'_standalone');
    if ~exist(modelFolder, 'dir')
        mkdir(modelFolder)
    end
    outputFile = append(modelFolder,'/',targSample,'.mat');
    model = SampleMod_Standalone.(targSample);
    save(outputFile,'model');
end

%% Get PAMat for reactions
PAMat_rxns = ones(numel(PanGEM.rxns),nSample);
for a = 1:nSample
    [junk, junk2, zeroOut_Idx] = intersect(canRemove{a},PanGEM.rxns);
    PAMat_rxns(zeroOut_Idx,a) = 0;
end

Meta_Assembly_L3 = Meta_Assembly;
Meta_Assembly_L3.PAMat_rxns = PAMat_rxns;

%% Save strain models and L3 structure
save(FileNames.Destination_L3,'Meta_Assembly_L3');

save(FileNames.Destination_SampleMod,'SampleMod');


%% OPTIONAL - Check strain models for growth under different limiting conditions
% load('data/models/StrMod.mat');
load(FileNames.Destination_SampleMod)
% assign limited models
limRxns = [{'AmmoniaEX'},{'OrthophosphateEX'},{'LightEX'}];
for a = 1:numel(limRxns)
    limRxns_idx(a) = find(strcmp(limRxns{a},PanGEM.rxns));
end
% check unlimited fluxes
unLimSol = solveLP(PanGEM,1);
unLim_limRxns_Flux = unLimSol.x(limRxns_idx);
% constrain limited model fluxes to half of optimal
nLim.lb(limRxns_idx(1)) = unLim_limRxns_Flux(1)./2;
pLim.lb(limRxns_idx(2)) = unLim_limRxns_Flux(2)./2;
LightLim.lb(limRxns_idx(3)) = unLim_limRxns_Flux(3)./2;


% check for growth
for a = 1:nSample;
    targSampleMod = SampleMod.(sample{a});
    % N lim
    b = 1;
    targSampleMod.lb(limRxns_idx(b)) = unLim_limRxns_Flux(b)./2;
    tempSol = solveLP(targSampleMod);
    if tempSol.stat
        NLimGrowth(a) = -tempSol.f;
    else
        NLimGrowth(a) = 0;
    end
    targSampleMod.lb(limRxns_idx(b)) = -1000; % reset to unbounded
    
    % P lim
    b = 2;
    targSampleMod.lb(limRxns_idx(b)) = unLim_limRxns_Flux(b)./2;
    tempSol = solveLP(targSampleMod);
    if tempSol.stat
        PLimGrowth(a) = -tempSol.f;
    else
        PLimGrowth(a) = 0;
    end
    targSampleMod.lb(limRxns_idx(b)) = -1000; % reset to unbounded
    
    % Light lim
    b = 3;
    targSampleMod.lb(limRxns_idx(b)) = unLim_limRxns_Flux(b)./2;
    tempSol = solveLP(targSampleMod);
    if tempSol.stat
        LightLimGrowth(a) = -tempSol.f;
    else
        LightLimGrowth(a) = 0;
    end
    targSampleMod.lb(limRxns_idx(b)) = -1000; % reset to unbounded
end









