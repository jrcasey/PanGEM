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

% (Less than 2 minutes runtime)
% Tatsuro Tanioka 20220725- Modified "Step_6_Strain_Models.m" by John Casey

%% 
clear

% Set transect here
% Transect = "C13"
% Transect = "I09"
% Transect = "I07"
% Transect = "Globe_AMT_P18_etc"

Transects = ["C13","I09","I07","Globe_AMT_P18_etc"];
% Confidence level (alpha): 90, 95, or 99%
alphas = {'alpha90', 'alpha95', 'alpha99'};

N = length(Transects);
% Pre-allocate the array of structures
for i = 1:N;
    SampleMod_Merged(i).Data = struct(Transects(i),[]);
end

%% Paths
working_dir  = '/Users/tatsurotanioka/Desktop/Project/BioGOSHIP_MetaPanGEM/';
FileNames.PanGEM_path = append(working_dir, 'PanGEM/data/models/Pro_PanGEM.mat');
FileNames.Destination_SampleMod_Merged_alpha90 = append(working_dir, 'Metagenome_GEM/models/SampleMod_Merged_alpha90.mat');
FileNames.Destination_SampleMod_Merged_alpha95 = append(working_dir, 'Metagenome_GEM/models/SampleMod_Merged_alpha95.mat');
FileNames.Destination_SampleMod_Merged_alpha99 = append(working_dir, 'Metagenome_GEM/models/SampleMod_Merged_alpha99.mat');
FileNames.Destination_Samplelist_Merged = append(working_dir, 'Metagenome_GEM/models/sampleList.mat');
%% Loop to Create SampleMod for each transect
for i = 1:length(Transects);
    clear SampleMod
    Transect = Transects(i)
    FileNames.L2Assembly_path = append(working_dir, 'Make_KO_Table/KO_PresAbs_Table/MATfiles/', Transect, '_Assembly_Uncurated.mat');
    FileNames.Destination_L3 = append(working_dir, 'Metagenome_GEM/Assembly_L3/',Transect,'_Assembly_L3','.mat');
    FileNames.Destination_SampleMod = append(working_dir, 'Metagenome_GEM/models/',Transect,'_SampleMod','.mat');
    FileNames.Destination_LimGrowth = append(working_dir, 'Metagenome_GEM/models/',Transect,'_LimGrowth','.mat');
    %% add PanGEM codes to the path
    addpath(genpath(append(working_dir,'PanGEM/src/')))
    %% Import PanGEM and the L2 Assembly
    load(FileNames.PanGEM_path)

    load(FileNames.L2Assembly_path)
    %% Use Uncurated Presence/Abs Matrix (i.e. CoreKOs are not forced and includes all KOs including those outside PanGEM KOs)

    Meta_Assembly = Meta_Assembly_Uncurated;
    sample = Meta_Assembly.SampleID_SCCGselected;
    % Remove deep samples from NH1418
    if Transect == "Globe_AMT_P18_etc";
        deep_samples = {'NH1418_134', 'NH1418_232', 'NH1418_328'};
        sample = setdiff(sample, deep_samples);    
    end
    nSample = numel(sample);

    %% Get essential reaction subset from PanGEM
    % Precalculated and saves as 'PanGEM/data/models/essentialRxns.mat'
    load(append(working_dir,'PanGEM/data/models/essentialRxns.mat'));
    % % BIOMASSMinimal contains the minimum set of BOF compounds required by all
    % % strains to produce biomass. Any additional components are considered
    % % non-essential
    % 
    % % minBOFind = find(strcmp('BIOMASSMinimal',PanGEM.rxns));
    % 
    % % BIOMASSCRUDE contains the full set of BOF compounds required by all
    % % strains to produce biomass. This will yield the most conservative
    % % essential subset.
    % minBOFind = find(strcmp('BIOMASSCRUDE',PanGEM.rxns));
    % 
    % PanGEM.c = zeros(numel(PanGEM.rxns),1);
    % PanGEM.c(minBOFind) = 1;
    % 
    % [essentialRxns, results] = solveCS(PanGEM);
    % % also don't mess with photosynthesis, oxphos and gap-filled reactions...
    % % so let's add those to the list
    % gapFilled = find(cellfun(@isempty,PanGEM.grRules));
    % PSInd = find(strcmp('Photosynthesis',PanGEM.subSystems));
    % OxPhosInd = find(strcmp('Oxidative phosphorylation',PanGEM.subSystems));
    % essentialRxns2 = unique([essentialRxns, PSInd',OxPhosInd',gapFilled']);
    % essentialRxns = essentialRxns2;
    % 
    % % specific (e.g., experimentally known) reactions to remove from this list
    % rmList = [{'R00475'}]; % glycolate oxidase
    % for a = 1:numel(rmList)
    %     rmIdx(a) = find(strcmp(rmList,PanGEM.rxns));
    % end
    % essentialRxns(find(ismember(essentialRxns,rmIdx))) = [];
    % 
    % save(append(working_dir,'PanGEM/data/models/essentialRxns.mat'),'essentialRxns');
    %% Get reactions to remove for each metagenome sample and for each confidence level
    nAlpha = numel(alphas);
    for a = 1:nSample
        for b = 1:nAlpha
            [rxnsToRemove{a,b}] = getRxnsToRemove_metagenomes(PanGEM,Meta_Assembly, sample{a}, alphas(b));
        end
    end

    %% Exclude essential reactions from reactions to remove and generate metagenome based models
    load(append(working_dir,'PanGEM/data/models/essentialRxns.mat'));

    for a = 1:nSample
        for b = 1:nAlpha
            tempRxnstoRemove = rxnsToRemove{a,b};
            canRemove{a,b} = setdiff(tempRxnstoRemove,PanGEM.rxns(essentialRxns));
        end
    end
    PAMat_model = zeros(numel(PanGEM.genes),nSample,nAlpha); 
    for a = 1:nSample
        for b = 1:nAlpha
            targSample = sample{a};
            [junk, removeIdx] = intersect(PanGEM.rxns,canRemove{a,b});
            SampleMod(b).(targSample) = PanGEM;
            % zero out S matrix columns corresponding to these reactions
            SampleMod(b).(targSample).S(:,removeIdx) = 0;
            % store geneset in a special model field (this is to avoid indexing
            % issues with actually removing genes)
            sampleRxns_Idx = find(sum(abs(SampleMod(b).(targSample).S),1));
            [junk, sampleGenes_Idx] = find(PanGEM.rxnGeneMat(sampleRxns_Idx,:));
            SampleMod(b).(targSample).geneSubset = unique(PanGEM.genes(sampleGenes_Idx));

            tempSol = solveLP(SampleMod(b).(targSample));
            if tempSol.stat
                mu1(a,b) = -tempSol.f;
            end
            [junk, include_idx] = intersect(PanGEM.genes,SampleMod(b).(targSample).geneSubset);
            PAMat_model(include_idx,a,b) = 1;
            SampleMod(b).PresenceAbsenceMatrix = PAMat_model(:,:,b);
        end
    end

    %% 
    modelFolder = append(working_dir, 'Metagenome_GEM/models/','SampleMod_standalone');
    if ~exist(modelFolder, 'dir')
        mkdir(modelFolder)
    end

    for b = 1:nAlpha
        for a = 1:nSample
            targSample = sample{a};
            % Get indices of genes and reactions to remove. 
            [junk, removeGenes_Idx] = setdiff(PanGEM.genes, SampleMod(b).(targSample).geneSubset);
            [junk, removeRxn_Idx] = intersect(PanGEM.rxns,canRemove{a,b});

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
            SampleMod_Standalone.(targSample).modelID = strcat('iTT',mat2str(numel(SampleMod(b).(targSample).geneSubset)),'_',targSample,'_',alphas{b});
            SampleMod_Standalone.(targSample).id = strcat(targSample,'_',alphas{b});

           % Export standalone model to mat
            outputFile = append(modelFolder,'/',targSample,'_',alphas{b},'.mat');
            model = SampleMod_Standalone.(targSample);
            save(outputFile,'model');
        end
    end

    %% Get PAMat for reactions
    PAMat_rxns = ones(numel(PanGEM.rxns),nSample,nAlpha);
    for a = 1:nSample
        for b = 1:nAlpha
            [junk, junk2, zeroOut_Idx] = intersect(canRemove{a,b},PanGEM.rxns);
            PAMat_rxns(zeroOut_Idx,a,b) = 0;
        end
    end

    Meta_Assembly_L3 = Meta_Assembly;
    Meta_Assembly_L3.PAMat_rxns_alpha90 = PAMat_rxns(:,:,1);
    Meta_Assembly_L3.PAMat_rxns_alpha95 = PAMat_rxns(:,:,2);
    Meta_Assembly_L3.PAMat_rxns_alpha99 = PAMat_rxns(:,:,3);
    %% Save metagenome sample models and L3 structure
    save(FileNames.Destination_L3,'Meta_Assembly_L3');

    SampleMod = rmfield(SampleMod, 'PresenceAbsenceMatrix');
    save(FileNames.Destination_SampleMod,'SampleMod');
    SampleMod_Merged(i).Data = SampleMod;
end
%% Merge all the Sample Mod from different cruises into one
SampleMod_C13 = SampleMod_Merged(1).Data;
SampleMod_I09 = SampleMod_Merged(2).Data;
SampleMod_I07 = SampleMod_Merged(3).Data;
SampleMod_Globe_AMT_P18_etc = SampleMod_Merged(4).Data;

mergestructs_4 = @(w,x,y,z) cell2struct([struct2cell(w);struct2cell(x);struct2cell(y);struct2cell(z)],[fieldnames(w);fieldnames(x);fieldnames(y);fieldnames(z)]);
SampleMod_Merged = mergestructs_4(SampleMod_I09,SampleMod_Globe_AMT_P18_etc,SampleMod_I07,SampleMod_C13);
sampleList = fieldnames(SampleMod_Merged);
save(FileNames.Destination_Samplelist_Merged,'sampleList');

% Save SampleMod Merged files separately for each alpha as they are very large
SampleMod_Merged_alpha90 = SampleMod_Merged(1);
SampleMod_Merged_alpha95 = SampleMod_Merged(2);
SampleMod_Merged_alpha99 = SampleMod_Merged(3);
save(FileNames.Destination_SampleMod_Merged_alpha90,'SampleMod_Merged_alpha90');
save(FileNames.Destination_SampleMod_Merged_alpha95,'SampleMod_Merged_alpha95');
save(FileNames.Destination_SampleMod_Merged_alpha99,'SampleMod_Merged_alpha99');
%% OPTIONAL - Check metagenome models for growth under different limiting conditions
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
for a = 1:nSample
    for b2 = 1:nAlpha
        targSampleMod = SampleMod(b2).(sample{a});
        % N lim
        b = 1;
        targSampleMod.lb(limRxns_idx(b)) = unLim_limRxns_Flux(b)./2;
        tempSol = solveLP(targSampleMod);
        if tempSol.stat
            NLimGrowth(a,b2) = -tempSol.f;
        else
            NLimGrowth(a,b2) = 0;
        end
        targSampleMod.lb(limRxns_idx(b)) = -1000; % reset to unbounded

        % P lim
        b = 2;
        targSampleMod.lb(limRxns_idx(b)) = unLim_limRxns_Flux(b)./2;
        tempSol = solveLP(targSampleMod);
        if tempSol.stat
            PLimGrowth(a,b2) = -tempSol.f;
        else
            PLimGrowth(a,b2) = 0;
        end
        targSampleMod.lb(limRxns_idx(b)) = -1000; % reset to unbounded

        % Light lim
        b = 3;
        targSampleMod.lb(limRxns_idx(b)) = unLim_limRxns_Flux(b)./2;
        tempSol = solveLP(targSampleMod);
        if tempSol.stat
            LightLimGrowth(a,b2) = -tempSol.f;
        else
            LightLimGrowth(a,b2) = 0;
        end
        targSampleMod.lb(limRxns_idx(b)) = -1000; % reset to unbounded
    end
end
% Save Growth Rates
LimGrowth.NLimGrowth = NLimGrowth;
LimGrowth.PLimGrowth = PLimGrowth;
LimGrowth.LightLimGrowth = LightLimGrowth;
save(FileNames.Destination_LimGrowth,'LimGrowth');





