%% Step 7 - Analysis
% Just a collection of a few useful plots and analyses to have on hand

% - Rarefaction curves for the pangenome/panreactome and the core subsets.
% - Gene/reaction/metabolite subsets for each metagenomes.
% - Growth potential on all available sole sources of C,N,P, and S.

% (About 20 minutes runtime)
% Tatsuro Tanioka for metganome-PanGEM(July 26, 2022)
% Original:"Step_7_Analysis.m" by John Casey

clear
close all
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
% alphas = {'alpha90', 'alpha95', 'alpha99'};
alphas = {'alpha90', 'alpha95', 'alpha99','top910','top920','top930','top940','top950'};

%% Matrices of KO gene presence/absence for each sample-derived GEMs

FileNames.L2Assembly_path = append(working_dir, 'Make_KO_Table/KO_PresAbs_Table/MATfiles/', Transect, '_Assembly.mat');
load(FileNames.L2Assembly_path)

FileNames.PAgenes_Mod = append(working_dir, 'Metagenome_GEM/models/',Transect,'_PAgenes_Mod','.mat');
load(FileNames.PAgenes_Mod);

% parse out sample names
samples = Meta_Assembly.SampleID;
if Transects{i} == "Globe_AMT_P18_etc";
    deep_samples = {'NH1418_134', 'NH1418_232', 'NH1418_328'};
    samples = setdiff(samples, deep_samples);    
end
nSample = numel(samples);
PAMat = PlotPresAbs_KO_GEM(PanGEM, PAgenes_Mod, samples, alphas, Transect, working_dir);


%% Matrices of KO gene presence/absence for each random-derived GEMs
totalKOs = [910, 920, 930, 940, 950];
totalKO_names = strcat("Random_top",string(totalKOs));
totalrndsamples = 100;
FileNames.Destination_PAgenes_RandomMod = append(working_dir, 'Metagenome_GEM/models/PAgenes_RandomMod_',num2str(totalrndsamples),'.mat');
load(FileNames.Destination_PAgenes_RandomMod)
rnd_names = strcat("rnd",string([1:totalrndsamples]));
PAMat_Rand_Mod = PlotPresAbs_KO_GEM(PanGEM, PAgenes_RandomMod,rnd_names, totalKO_names, "Random", working_dir);

%% Presence/Absence Table for all KOs recruited
% FileNames.Assembly_L3 = append(working_dir,"Make_KO_Table/KO_PresAbs_Table/MATfiles/",Transect,"_Assembly.mat");
% load(FileNames.Assembly_L3);
% 
% [PAMat] = PlotPresAbs_KO_all(Meta_Assembly, alphas, Transect, working_dir);

%% Function to visualize KO gene presence/absence for each sample-derived GEMs
function [PAMat] = PlotPresAbs_KO_GEM(PanGEM, PAgenes_Mod, samples, alphas, Transect, working_dir);
    nSample = numel(samples);
    nAlpha = numel(alphas);
    sample_nKOs = zeros(nSample,nAlpha);
    PAMat = zeros(numel(PanGEM.genes),nSample,nAlpha);
    for a = 1:nAlpha;
        alpha = alphas{a};
        PresenceAbsenceMatrix = PAgenes_Mod(a).PresenceAbsenceMatrix;        
        sample_nKOs(:,a) = sum(PresenceAbsenceMatrix,1);
        [KO_nSample] = sum(PresenceAbsenceMatrix,2);
        PAMat(:,:,a) = squeeze(PresenceAbsenceMatrix);
    end
    
    for a = 1:numel(alphas);
        f = figure;
        % ax2 = subplot(3,1,a);
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
        filename = append(working_dir,"Metagenome_GEM/Analysis_Plots/",Transect,"_",alphas(a),"_GEM_KO.pdf");
        orient(f, 'landscape')
        print(f, filename, '-dpdf', '-bestfit')
    end
   
end

%% Function to visualize all KO gene presence/absence from coverage table
% Red = Samples with SCCG coverage > 5
function [PAMat2] = PlotPresAbs_KO_all(Meta_Assembly, alphas, Transect, working_dir);

    samples = Meta_Assembly.SampleID;
    PresenceAbsenceMatrix(:,:,1) = Meta_Assembly.PresenceAbsenceMatrix_alpha90;
    PresenceAbsenceMatrix(:,:,2) = Meta_Assembly.PresenceAbsenceMatrix_alpha95;
    PresenceAbsenceMatrix(:,:,3) = Meta_Assembly.PresenceAbsenceMatrix_alpha99;    
    PresenceAbsenceMatrix(:,:,4) = Meta_Assembly.PresenceAbsenceMatrix_top910;
    PresenceAbsenceMatrix(:,:,5) = Meta_Assembly.PresenceAbsenceMatrix_top920;
    PresenceAbsenceMatrix(:,:,6) = Meta_Assembly.PresenceAbsenceMatrix_top930;
    PresenceAbsenceMatrix(:,:,7) = Meta_Assembly.PresenceAbsenceMatrix_top940;   
    PresenceAbsenceMatrix(:,:,8) = Meta_Assembly.PresenceAbsenceMatrix_top950;   
    nSample = numel(samples);
    KO_all =  Meta_Assembly.uniqueKO;
    CoreKO_PangeneCluster = readtable("~/Desktop/Project/BioGOSHIP_MetaPanGEM/Gene_match_PanGEM/Core_KO_PanGEM.txt");
    PAMat = PresenceAbsenceMatrix;
    
    f = figure;
    f.Position = [1200 800 633 1200];
    cmap = [1 1 1; 0 0.5 0; 1 0 0];
    colormap(cmap)
    for k = 1:3
        subplot(3,1,k)
        [str_nKOs, str_order] = sort(sum(PresenceAbsenceMatrix(:,:,k),1),'descend');
        [KO_nStr, KO_order] = sort(sum(PresenceAbsenceMatrix(:,:,k),2),'descend');
        % Pickout samples with SCCG coverage > 5
        HQ_strIdx = find(ismember(Meta_Assembly.SampleID,Meta_Assembly.SampleID));
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
        filename = append(working_dir,"Metagenome_GEM/Analysis_Plots/",Transect,"_alphas_KO_coverage.pdf");
        title(title_plot,'Interpreter', 'none')
        hold on
    end
    orient(f, 'portrait')
    print(f, filename, '-dpdf', '-bestfit')
    
    f = figure;
    f.Position = [1200 800 633 1200];
    cmap = [1 1 1; 0 0.5 0; 1 0 0];
    colormap(cmap)
    for k = 4:2:8
        i = (k-2)/2;
        subplot(3,1,i)
        [str_nKOs, str_order] = sort(sum(PresenceAbsenceMatrix(:,:,k),1),'descend');
        [KO_nStr, KO_order] = sort(sum(PresenceAbsenceMatrix(:,:,k),2),'descend');
        % Pickout samples with SCCG coverage > 5
        HQ_strIdx = find(ismember(Meta_Assembly.SampleID,Meta_Assembly.SampleID));
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
        filename = append(working_dir,"Metagenome_GEM/Analysis_Plots/",Transect,"_topKOs_KO_coverage.pdf");
        title(title_plot,'Interpreter', 'none')
        hold on
    end
    orient(f, 'portrait')
    print(f, filename, '-dpdf', '-bestfit')    
    
end