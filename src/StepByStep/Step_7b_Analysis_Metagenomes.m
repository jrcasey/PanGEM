%% Step 7 - Analysis
% Just a collection of a few useful plots and analyses to have on hand

% - Rarefaction curves for the pangenome/panreactome and the core subsets.
% - Gene/reaction/metabolite subsets for each metagenomes.
% - Growth potential on all available sole sources of C,N,P, and S.

% (About 20 minutes runtime)
% Tatsuro Tanioka for metganome-PanGEM(July 26, 2022)
% Original:"Step_7_Analysis.m" by John Casey

clear
%% Add PanGEM to the path
addpath(genpath("/Users/tatsurotanioka/Desktop/Project/BioGOSHIP_MetaPanGEM/PanGEM"));
%% Directories and load PanGEM
working_dir  = '/Users/tatsurotanioka/Desktop/Project/BioGOSHIP_MetaPanGEM/';

FileNames.PanGEM = append(working_dir, 'PanGEM/data/models/Pro_PanGEM.mat');
FileNames.SoleElements = append(working_dir, 'PanGEM/data/dbLinks/SoleElementSources.csv');
load(FileNames.PanGEM);

%% Set transect here
% 
Transect = "C13"
% Transect = "I09"
% Transect = "I07"
% Transect = "Globe_AMT_P18_etc"

Transects = {'C13','I09','I07','Globe_AMT_P18_etc'};
nTransects = numel(Transects);
% nTransects = 1;
alphas = {'alpha90', 'alpha95', 'alpha99'};


%% Matrices of KO gene presence/absence for each sample-derived GEMs
for i = 1:nTransects;

    FileNames.SampleMod = append(working_dir, 'Metagenome_GEM/models/',Transects{i},'_SampleMod','.mat');
    load(FileNames.SampleMod);

    % parse out sample names
    samples = fieldnames(SampleMod);
    samples(ismember(samples,'PresenceAbsenceMatrix')) = [];
    nSample = numel(samples);

    PAMat = PlotPresAbs_KO_GEM(PanGEM, SampleMod, samples, alphas, Transects{i}, working_dir);
end
%% Presence/Absence Table for all KOs recruited
FileNames.Assembly_L3 = append(working_dir,"Make_KO_Table/KO_PresAbs_Table/MATfiles/",Transect,"_Assembly_Uncurated.mat");
load(FileNames.Assembly_L3);

[PAMat] = PlotPresAbs_KO_all(Meta_Assembly_Uncurated, alphas, Transect, working_dir);
%% Check growth potential for all strains
% This last section takes the longest time.

% list of elements to test
elements = [{'C'},{'N'},{'P'},{'S'}];
%elements = [{'N'},{'P'},{'S'}];

% import list of exomets and their corresponding uptake reaction IDs 
soleElements = readtable(FileNames.SoleElements,'Delimiter',',','ReadVariableNames',true);

% generate structure with fields for each element
soleSources = struct;

% get metabolites for each element
for a = 1:numel(elements)
    element_Idx = find(strcmp(elements{a},soleElements.Element));
    soleSources.(elements{a}).element_Idx = element_Idx;
    soleSources.(elements{a}).mets = soleElements.Met(element_Idx);
    soleSources.(elements{a}).TpRxn = soleElements.Uptake_rxn(element_Idx);
    soleSources.(elements{a}).Growth = zeros(numel(element_Idx),nStr);
end

% Generate models deficient in each element, and check that they are infeasible. 

% For carbon, prevent CO2 fixation uptake
RubiscoTransIdx = find(strcmp('R00024',PanGEM.rxns));

% For nitrogen, prevent ammonia uptake
AmmoniaTransIdx = find(strcmp('AmmoniaTRANS',PanGEM.rxns));

% for phosphorus, prevent orthophosphate uptake
PiTransIdx = find(strcmp('OrthophosphateTRANS',PanGEM.rxns));

% for sulfur, prevent sulfate uptake
SO4TransIdx = find(strcmp('SulfateTRANS',PanGEM.rxns));

% now loop through each strain and check for growth on each element source
threshold = 1e-3;
for a = 1:nStr
    % initial model
    tempMod = StrMod.(strains{a});
    for b = 1:numel(elements)
        % assign element limited model, loop through candidate sources, and
        % store output
        limMod = tempMod;
        if elements{b} == 'C'
            limMod.ub(RubiscoTransIdx) = 0;
            nExoMets = numel(soleSources.(elements{b}).mets);
            for c = 1:nExoMets
                limMod2 = limMod;
                TpRxn_idx = find(strcmp(soleSources.(elements{b}).TpRxn{c},limMod.rxns));
                limMod2.lb(TpRxn_idx) = -1000;
                tempSol = solveLP(limMod2);
                if tempSol.stat
                    if tempSol.f < -threshold
                        soleSources.(elements{b}).Growth(c,a) = 1;
                    end
                end
            end
        end
        
        if elements{b} == 'N'
            limMod.lb(AmmoniaTransIdx) = 0;
            nExoMets = numel(soleSources.(elements{b}).mets);
            for c = 1:nExoMets
                limMod2 = limMod;
                TpRxn_idx = find(strcmp(soleSources.(elements{b}).TpRxn{c},limMod.rxns));
                limMod2.lb(TpRxn_idx) = -1000;
                tempSol = solveLP(limMod2);
                if tempSol.stat
                    if tempSol.f < -threshold
                        soleSources.(elements{b}).Growth(c,a) = 1;
                    end
                end
            end
        end
        if elements{b} == 'P'
            limMod.lb(PiTransIdx) = 0;
            nExoMets = numel(soleSources.(elements{b}).mets);
            for c = 1:nExoMets
                limMod2 = limMod;
                TpRxn_idx = find(strcmp(soleSources.(elements{b}).TpRxn{c},limMod.rxns));
                limMod2.lb(TpRxn_idx) = -1000;
                tempSol = solveLP(limMod2);
                if tempSol.stat
                    if tempSol.f < -threshold
                        soleSources.(elements{b}).Growth(c,a) = 1;
                    end
                end
            end
        end
        if elements{b} == 'S'
            limMod.lb(SO4TransIdx) = 0;
            nExoMets = numel(soleSources.(elements{b}).mets);
            for c = 1:nExoMets
                limMod2 = limMod;
                TpRxn_idx = find(strcmp(soleSources.(elements{b}).TpRxn{c},limMod.rxns));
                limMod2.lb(TpRxn_idx) = -1000;
                tempSol = solveLP(limMod2);
                if tempSol.stat
                    if tempSol.f < -threshold
                        soleSources.(elements{b}).Growth(c,a) = 1;
                    end
                end
            end
        end
    end
end


% plot results
% Group strains by ecotype
ecotypes = [{'HLI'},{'HLII'},{'LLI'},{'LLII_LLIII'},{'LLIV'}];
for a = 1:nStr
    tempIdx = find(strcmp(strains{a},Pro_Assembly_L3.orgDatabase.StrainName));
    StrIdx(a) = tempIdx(1);
end
for a = 1:numel(ecotypes)
    ecotype_idx{a} = find(strcmp(ecotypes{a},Pro_Assembly_L3.orgDatabase.Ecotype(StrIdx)));
end
ecotype_ordered = vertcat(ecotype_idx{:});

allSources_Growth = [soleSources.C.Growth; soleSources.N.Growth; soleSources.P.Growth; soleSources.S.Growth];
%allSources_Growth = [soleSources.N.Growth; soleSources.P.Growth; soleSources.S.Growth];

figure
imagesc(allSources_Growth(:,ecotype_ordered));
set(gca,'XTick',1:nStr,'XTickLabels',strains(ecotype_ordered));
%set(gca,'YTick',1:numel(soleElements.Met([soleSources.N.element_Idx; soleSources.P.element_Idx; soleSources.S.element_Idx])),'YTickLabels',strrep(soleElements.Met([soleSources.N.element_Idx; soleSources.P.element_Idx; soleSources.S.element_Idx]),'_','-'))
set(gca,'YTick',1:numel(soleElements.Met),'YTickLabels',strrep(soleElements.Met,'_','-'))
set(gca,'FontSize',10)
ylabel('Sole source')
xlabel('Strain')
xtickangle(90)


%% Function to visualize KO gene presence/absence for each sample-derived GEMs
function [PAMat] = PlotPresAbs_KO_GEM(PanGEM, SampleMod, samples, alphas, Transect, working_dir);
    nSample = numel(samples);
    nAlpha = numel(alphas);
    sample_nKOs = zeros(nSample,nAlpha);
    PAMat = zeros(numel(PanGEM.genes),nSample,nAlpha);
    f = figure;
    f.Position = [1200 800 633 1200];
    for a = 1:nAlpha;
        alpha = alphas{a};
        PresenceAbsenceMatrix = SampleMod(a).PresenceAbsenceMatrix;        
        sample_nKOs(:,a) = sum(PresenceAbsenceMatrix,1);
        [KO_nSample] = sum(PresenceAbsenceMatrix,2);
        PAMat(:,:,a) = squeeze(PresenceAbsenceMatrix);
    end

    for a = 1:nAlpha;
        ax2 = subplot(3,1,a);
        imagesc(PAMat(:,:,a));
        colormap(flipud(hot));
        %title('All Strains')
        xticks(1:1:length(samples));
        xtickangle(90)
        xlabel_name = samples;
        xticklabels(xlabel_name(1:1:end));
        xlim([1 width(PAMat)])
        % xlabel('Sample')
        ylabel('GEM_KO','Interpreter', 'none')
        set(gca,'FontSize',7)
        set(gca,'TickLabelInterpreter','none')
        title_plot=append(Transect,"_GEM_KO_",alphas{a});
        title(title_plot,'Interpreter', 'none')
        h=gca; h.YAxis.TickLength = [0 0];
        h=gca; h.XAxis.TickLength = [0 0];
        h2 = text(1:1:length(samples), repelem(length(KO_nSample)+350,length(1:1:length(samples))),string(sample_nKOs(:,a)), "FontSize",6,'Interpreter', 'none','FontWeight','bold');
        set(h2,'Rotation',90);
        text(-10, (length(KO_nSample))+300, "Total KO", "FontSize",7,'Interpreter', 'none')
    hold on
    end
    filename = append(working_dir,"Metagenome_GEM/Analysis_Plots/",Transect,"_GEM_KO.pdf");
    orient(f, 'portrait')
    print(f, filename, '-dpdf', '-bestfit')

end

%% Function to visualize all KO gene presence/absence from coverage table
% Red = Samples with SCCG coverage > 5
function [PAMat2] = PlotPresAbs_KO_all(Meta_Assembly_Uncurated, alphas, Transect, working_dir);

    samples = Meta_Assembly_Uncurated.SampleID;
    PresenceAbsenceMatrix(:,:,1) = Meta_Assembly_Uncurated.PresenceAbsenceMatrix_alpha90;
    PresenceAbsenceMatrix(:,:,2) = Meta_Assembly_Uncurated.PresenceAbsenceMatrix_alpha95;
    PresenceAbsenceMatrix(:,:,3) = Meta_Assembly_Uncurated.PresenceAbsenceMatrix_alpha99;    
    nSample = numel(samples);
    KO_all =  Meta_Assembly_Uncurated.uniqueKO;
    CoreKO_PangeneCluster = readtable("~/Desktop/Project/BioGOSHIP_MetaPanGEM/Gene_match_PanGEM/Core_KO_PanGEM.txt");
    PAMat = PresenceAbsenceMatrix;
    f = figure;
    f.Position = [1200 800 633 1200];
    cmap = [1 1 1; 0 0.5 0; 1 0 0];
    colormap(cmap)
    for k = 1:numel(alphas)
        subplot(3,1,k)
        [str_nKOs, str_order] = sort(sum(PresenceAbsenceMatrix(:,:,k),1),'descend');
        [KO_nStr, KO_order] = sort(sum(PresenceAbsenceMatrix(:,:,k),2),'descend');
        % Pickout samples with SCCG coverage > 5
        HQ_strIdx = find(ismember(Meta_Assembly_Uncurated.SampleID,Meta_Assembly_Uncurated.SampleID_SCCGselected));
        for a = 1:numel(HQ_strIdx)
            HQ_strIdx2(a) = find(HQ_strIdx(a) == str_order);
        end
        PAMat2 = PresenceAbsenceMatrix(KO_order,str_order,k);
        PAMat2(:,HQ_strIdx2) = PAMat2(:,HQ_strIdx2)*2;
        samples2 = samples(str_order);
        imagesc(PAMat2);
        xlim([1 nSample])
        ylim([0 5000]);
        % xlabel('Sample')
        xlabel_name = samples2;
        xticks(1:4:length(samples2));
        xticklabels(xlabel_name(1:3:end));
        yticks(0:500:5000);
        set(gca,'TickLabelInterpreter','none')
        ylabel('KO')
        xtickangle(90)
        h=gca; h.YAxis.TickLength = [0 0];
        h=gca; h.XAxis.TickLength = [0 0];
        set(gca,'FontSize',8)
        h2 = text(1:4:length(samples2), repelem(length(KO_nStr)+1400,length(1:4:length(samples2))),string(str_nKOs(1:4:end)), "FontSize",6,'Interpreter', 'none','FontWeight','bold');
        set(h2,'Rotation',90);
        text(-20, (length(KO_nStr))+1100, "Total KO", "FontSize",7,'Interpreter', 'none')
        title_plot=append(Transect,"_KO_coverage_",alphas{k});
        filename = append(working_dir,"Metagenome_GEM/Analysis_Plots/",Transect,"_KO_coverage.pdf");
        title(title_plot,'Interpreter', 'none')
        hold on
    end
    orient(f, 'portrait')
    print(f, filename, '-dpdf', '-bestfit')
end