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

%% 
% Set transect here
Transect = "C13"
% Transect = "I09"
% Transect = "I07"
% Transect = "Globe_AMT_P18_etc"
%% Paths
working_dir  = '/Users/tatsurotanioka/Desktop/Project/BioGOSHIP_MetaPanGEM/';
FileNames.PanGEM_path = append(working_dir, 'PanGEM/data/models/Pro_PanGEM.mat');
FileNames.L2Assembly_path = append(working_dir, 'Make_KO_Table/KO_PresAbs_Table/MATfiles/', Transect, '_Assembly.mat');
FileNames.Destination_L3 = append(working_dir, 'Metagenome_GEM/',Transect,'_Assembly_L3.mat');
FileNames.Destination_StrMod = append(working_dir, 'Metagenome_GEM/models/',Transect,'_StrMod.mat');
FileNames.Destination_StrMod_Standalone = append(working_dir, 'Metagenome_GEM/models/',Transect,'_StrMod_Standalone.mat');

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

%% Get essential reaction subset
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
essentialRxns2 = unique([essentialRxns, PSInd',OxPhosInd',gapFilled'])
essentialRxns = essentialRxns2;

% specific (e.g., experimentally known) reactions to remove from this list
rmList = [{'R00475'}]; % glycolate oxidase
for a = 1:numel(rmList)
    rmIdx(a) = find(strcmp(rmList,PanGEM.rxns));
end
essentialRxns(find(ismember(essentialRxns,rmIdx))) = [];

%% Get reactions to remove for each metagenome sample
   
for a = 1:nStrains
    [rxnsToRemove{a}] = getRxnsToRemove(PanGEM,Pro_Assembly_L2, strain{a});
end
%%
for a = 1:nSample
    [rxnsToRemove{a}] = getRxnsToRemove(PanGEM,Meta_Assembly, sample{a});
end

%% Exclude essential reactions from reactions to remove and generate strain models
for a = 1:nStrains
    tempRxnstoRemove = rxnsToRemove{a};
    canRemove{a} = setdiff(tempRxnstoRemove,PanGEM.rxns(essentialRxns));
end

for a = 1:nStrains
    targStr = strain{a};
    [junk, removeIdx] = intersect(PanGEM.rxns,canRemove{a});
    StrMod.(targStr) = PanGEM;
    % zero out S matrix columns corresponding to these reactions
    StrMod.(targStr).S(:,removeIdx) = 0;
    % store geneset in a special model field (this is to avoid indexing
    % issues with actually removing genes)
    strRxns_Idx = find(sum(abs(StrMod.(targStr).S),1));
    [junk, strGenes_Idx] = find(PanGEM.rxnGeneMat(strRxns_Idx,:));
    StrMod.(targStr).geneSubset = unique(PanGEM.genes(strGenes_Idx));
    
    tempSol = solveLP(StrMod.(targStr));
    if tempSol.stat
        mu1(a) = -tempSol.f;
    end
end

%% Generate standalone strain models
for a = 1:nStrains
    targStr = strain{a};
    % Get indices of genes and reactions to remove. 
    [junk, removeGenes_Idx] = setdiff(PanGEM.genes, StrMod.(targStr).geneSubset);
    [junk, removeRxn_Idx] = intersect(PanGEM.rxns,canRemove{a});
    
    % Initialize strain model
    StrMod_Standalone.(targStr) = PanGEM;
    
    % Gene fields
    StrMod_Standalone.(targStr).genes(removeGenes_Idx) = [];
    StrMod_Standalone.(targStr).rxnGeneMat(:,removeGenes_Idx) = [];
    StrMod_Standalone.(targStr).geneMiriams(removeGenes_Idx) = [];
    StrMod_Standalone.(targStr).geneSBOTerms(removeGenes_Idx) = [];
    
    % Reaction fields
    StrMod_Standalone.(targStr).S(:,removeRxn_Idx) = [];
    StrMod_Standalone.(targStr).rxns(removeRxn_Idx) = [];
    StrMod_Standalone.(targStr).lb(removeRxn_Idx) = [];
    StrMod_Standalone.(targStr).ub(removeRxn_Idx) = [];
    StrMod_Standalone.(targStr).c(removeRxn_Idx) = [];
    StrMod_Standalone.(targStr).rules(removeRxn_Idx) = [];
    StrMod_Standalone.(targStr).grRules(removeRxn_Idx) = [];
    StrMod_Standalone.(targStr).rxnGeneMat(removeRxn_Idx,:) = [];
    StrMod_Standalone.(targStr).rxnNames(removeRxn_Idx) = [];
    StrMod_Standalone.(targStr).subSystems(removeRxn_Idx) = [];
    StrMod_Standalone.(targStr).eccodes(removeRxn_Idx) = [];
    StrMod_Standalone.(targStr).rxnComps(removeRxn_Idx) = [];
    StrMod_Standalone.(targStr).rxnKEGGID(removeRxn_Idx) = [];
    StrMod_Standalone.(targStr).rxnMetaNetXID(removeRxn_Idx) = [];
    StrMod_Standalone.(targStr).rxnSBOTerms(removeRxn_Idx) = [];
    StrMod_Standalone.(targStr).rev(removeRxn_Idx) = [];
    
    StrMod_Standalone.(targStr).rxnKEGGID = StrMod_Standalone.(targStr).rxnKEGGID';
    % Remove unused metabolites
    unusedMet_Idx = find(sum(abs(StrMod_Standalone.(targStr).S),2)==0);

    % Metabolite fields
    StrMod_Standalone.(targStr).S(unusedMet_Idx,:) = [];
    StrMod_Standalone.(targStr).mets(unusedMet_Idx) = [];
    StrMod_Standalone.(targStr).b(unusedMet_Idx) = [];
    StrMod_Standalone.(targStr).csense(unusedMet_Idx) = [];
    StrMod_Standalone.(targStr).metCharges(unusedMet_Idx) = [];
    StrMod_Standalone.(targStr).metFormulas(unusedMet_Idx) = [];
    StrMod_Standalone.(targStr).metNames(unusedMet_Idx) = [];
    StrMod_Standalone.(targStr).metComps(unusedMet_Idx) = [];    
    StrMod_Standalone.(targStr).metKEGGID(unusedMet_Idx) = [];
    StrMod_Standalone.(targStr).metPubChemID(unusedMet_Idx) = [];
    StrMod_Standalone.(targStr).metChEBIID(unusedMet_Idx) = [];
    StrMod_Standalone.(targStr).metInChI(unusedMet_Idx) = [];
    StrMod_Standalone.(targStr).metMetaNetXID(unusedMet_Idx) = [];
    StrMod_Standalone.(targStr).metSBOTerms(unusedMet_Idx) = [];
    
    % Annotate model field
    StrMod_Standalone.(targStr).description = targStr;
    StrMod_Standalone.(targStr).modelID = strcat('iJC',mat2str(numel(StrMod.(targStr).geneSubset)),'_',targStr);
    StrMod_Standalone.(targStr).id = targStr;

    % Export standalone model to mat
    outputFile = strcat('data/models/strains/',targStr,'.mat');
    model = StrMod_Standalone.(targStr);
    save(outputFile,'model');
end

%% Get PAMat for reactions
PAMat_rxns = ones(numel(PanGEM.rxns),nStrains);
for a = 1:nStrains
    [junk, junk2, zeroOut_Idx] = intersect(canRemove{a},PanGEM.rxns);
    PAMat_rxns(zeroOut_Idx,a) = 0;
end

Pro_Assembly_L3 = Pro_Assembly_L2;
Pro_Assembly_L3.PAMat_rxns = PAMat_rxns;

%% Save strain models and L3 structure
save(FileNames.Destination_L3,'Pro_Assembly_L3');

save(FileNames.Destination_StrMod,'StrMod');


%% OPTIONAL - Check strain models for growth under different limiting conditions
load('data/models/StrMod.mat');

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
for a = 1:nStrains
    targStrMod = StrMod.(strain{a});
    % N lim
    b = 1;
    targStrMod.lb(limRxns_idx(b)) = unLim_limRxns_Flux(b)./2;
    tempSol = solveLP(targStrMod);
    if tempSol.stat
        NLimGrowth(a) = -tempSol.f;
    else
        NLimGrowth(a) = 0;
    end
    targStrMod.lb(limRxns_idx(b)) = -1000; % reset to unbounded
    
    % P lim
    b = 2;
    targStrMod.lb(limRxns_idx(b)) = unLim_limRxns_Flux(b)./2;
    tempSol = solveLP(targStrMod);
    if tempSol.stat
        PLimGrowth(a) = -tempSol.f;
    else
        PLimGrowth(a) = 0;
    end
    targStrMod.lb(limRxns_idx(b)) = -1000; % reset to unbounded
    
    % Light lim
    b = 3;
    targStrMod.lb(limRxns_idx(b)) = unLim_limRxns_Flux(b)./2;
    tempSol = solveLP(targStrMod);
    if tempSol.stat
        LightLimGrowth(a) = -tempSol.f;
    else
        LightLimGrowth(a) = 0;
    end
    targStrMod.lb(limRxns_idx(b)) = -1000; % reset to unbounded
end









