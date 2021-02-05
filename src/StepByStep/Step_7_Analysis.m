%% Step 7 - Analysis
% Just a collection of a few useful plots and analyses to have on hand

% - Rarefaction curves for the pangenome/panreactome and the core subsets.
% - Gene/reaction/metabolite subsets for each clade of organisms.
% - Growth potential on all available sole sources of C,N,P, and S.

% (About 20 minutes runtime)

%% Directories
FileNames.Assembly_L3 = 'data/assemblies/Assembly_L3_20200830.mat';
FileNames.PanGEM = 'data/models/Pro_PanGEM.mat';
FileNames.StrMod = 'data/models/StrMod.mat';
FileNames.SoleElements = 'data/dbLinks/SoleElementSources.csv';

%% Load L3 Assembly, PanGEM, and StrMod
load(FileNames.Assembly_L3);
load(FileNames.PanGEM);
load(FileNames.StrMod);

% parse out strain names
strains = fieldnames(StrMod);
nStr = numel(strains);

%% Matrices of gene, reaction, and metabolite presence/absence
% All strains 
% rank both dimensions
[str_nKOs, str_order] = sort(sum(Pro_Assembly_L3.PresenceAbsenceMatrix,1),'descend');
[KO_nStr, KO_order] = sort(sum(Pro_Assembly_L3.PresenceAbsenceMatrix,2),'descend');

HQ_strIdx = find(Pro_Assembly_L3.orgDatabase.Include);
for a = 1:numel(HQ_strIdx)
    HQ_strIdx2(a) = find(HQ_strIdx(a) == str_order);
end

PAMat2 = Pro_Assembly_L3.PresenceAbsenceMatrix(KO_order,str_order);
PAMat2(:,HQ_strIdx2) = PAMat2(:,HQ_strIdx2)*2;


cmap = [1 1 1; 0 0.5 0; 1 0 0];

figure
ax1 = subplot(2,1,1)
plot(Pro_Assembly_L3.orgDatabase.Size_Mb(str_order),'-','Color',cmap(2,:),'LineWidth',3)
hold on
plot(HQ_strIdx2,Pro_Assembly_L3.orgDatabase.Size_Mb(str_order(HQ_strIdx2)),'.','MarkerSize',10,'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor',cmap(3,:))
set(gca,'XTick',[])
xlim([1 649])
ylabel('Genome Size [Mb]')
set(gca,'FontSize',20)
ax2 = subplot(2,1,2)
imagesc(PAMat2);
colormap(cmap)
%title('All Strains')
xlim([1 649])
xlabel('Strain')
ylabel('KO')
set(gca,'FontSize',20)

% There is a possiblity that strains AG_418_M21 and B245a_518D8 are not actually Pro...




%% Rarefaction curves (KO's)

% issues with contamination on indices 2 and 286
Pro_Assembly_L3.PresenceAbsenceMatrix(:,2) = Pro_Assembly_L3.PresenceAbsenceMatrix(:,285);
Pro_Assembly_L3.PresenceAbsenceMatrix(:,286) = Pro_Assembly_L3.PresenceAbsenceMatrix(:,285);

% Pan
nStrTotal = numel(Pro_Assembly_L3.orgDatabase.StrainName);
for a = 1:100
    randIdx = randi(nStrTotal,nStrTotal,1);
    cumulativeKOs{1} = Pro_Assembly_L3.uniqueKO(find(Pro_Assembly_L3.PresenceAbsenceMatrix(:,randIdx(1))));
    for b = 1:nStrTotal
        if b > 1
            newKOs = Pro_Assembly_L3.uniqueKO(find(Pro_Assembly_L3.PresenceAbsenceMatrix(:,randIdx(b))));
            cumulativeKOs{b} = unique([cumulativeKOs{b-1}; newKOs]);
            nKOs(b,a) = numel(cumulativeKOs{b});
        end
    end
end

figure
h1 = plot(2:nStrTotal,mean(nKOs(2:end,:),2),'-r','LineWidth',3);
hold on
h2 = errorbar(2:nStrTotal,mean(nKOs(2:end,:),2),std(nKOs(2:end,:),[],2),'-k');
set(gca,'FontSize',20)
xlabel('Number of genomes')
ylabel('Number of KOs')

% Core
% indices of core strains
for a = 1:nStr
    tempIdx = find(strcmp(strains{a},Pro_Assembly_L3.orgDatabase.StrainName));
    strIdx(a) = tempIdx(1);
end

for a = 1:100
    randIdx = randi(nStr,nStr,1);
    cumulativeKOs{1} = Pro_Assembly_L3.uniqueKO(find(Pro_Assembly_L3.PresenceAbsenceMatrix(:,strIdx(randIdx(1)))));
    for b = 1:nStr
        if b > 1
            newKOs = Pro_Assembly_L3.uniqueKO(find(Pro_Assembly_L3.PresenceAbsenceMatrix(:,strIdx(randIdx(b)))));
            cumulativeKOs{b} = unique(intersect(cumulativeKOs{b-1},newKOs));
            nCoreKOs(b,a) = numel(cumulativeKOs{b});
        end
    end
end

figure
h1 = plot(2:nStr,mean(nCoreKOs(2:end,:),2),'-k','LineWidth',3);
hold on
h2 = errorbar(2:nStr,mean(nCoreKOs(2:end,:),2),std(nCoreKOs(2:end,:),[],2),'-k');
set(gca,'FontSize',20)
xlabel('Number of genomes')
ylabel('Number of KOs')

% both
figure
subplot(2,1,1)
h1 = plot(2:nStrTotal,mean(nKOs(2:end,:),2),'-r','LineWidth',3);
hold on
h2 = errorbar(2:nStrTotal,mean(nKOs(2:end,:),2),std(nKOs(2:end,:),[],2),'-k');
set(h2,'LineStyle','none')
ylim([700 2100])
set(gca,'FontSize',20)
xlabel('Number of genomes')
ylabel('Number of KOs')
subplot(2,1,2)
h3 = plot(2:nStr,mean(nCoreKOs(2:end,:),2),'-r','LineWidth',3);
hold on
h4 = errorbar(2:nStr,mean(nCoreKOs(2:end,:),2),std(nCoreKOs(2:end,:),[],2),'-k');
set(h4,'LineStyle','none')
ylim([750 950])
set(gca,'FontSize',20)
xlabel('Number of genomes')
ylabel('Number of KOs')




%% Check growth potential for all strains

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
    tempIdx = find(strcmp(strains{a},Pro_Assembly_L3.orgDatabase.StrainName))
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

