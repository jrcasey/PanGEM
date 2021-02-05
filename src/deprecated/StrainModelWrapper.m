%% Generate strain models
% Import model
fileName = sprintf('%s','CBIOMES/Pangenomes/Prochlorococcus/PanGEM_pro_20200105.xlsx');
PanGEM = importExcelModel(fileName,true,true,false);
[sol hsSolOut] = solveLP(PanGEM,1);
% optional: set objective to the minimum biomass composition BOF
% minBOFind = find(strcmp('BIOMASSMinimal',PanGEM.rxns));
% PanGEM.c = zeros(numel(PanGEM.rxns),1);
% PanGEM.c(minBOFind) = 1;
% sol2 = solveLP(PanGEM,1);

load('Pro_Assembly3.mat');
strain = Dat.orgDatabase.Strain(Dat.orgDatabase.Scaffolds<50);
nStrains = numel(strain);

%% Algorithm 2 (random search)
% This one takes for.ev.er
% nPermutations = 100;
% minGrowth = 0.1;
% scaffoldCutoff = 50;
% fileName = 'CBIOMES/Pangenomes/Prochlorococcus/AllStrainsMinBOF4.mat';
% [targStrMod nRxns2 targStrGrowth targStrFluxes targStrShadow] = getStrainMod_Algorithm2_v1(Dat,PanGEM,nPermutations,minGrowth,scaffoldCutoff,fileName);
% load AllStrainsMinBOF4.mat
% 
% 
% for i = 1:numel(strain)
%     nRxnsFinal(i) = numel(targStrMod.(strain{i}).rxns);
% end


%% L1-norm minimization (Compressive sensing) algorithm
% optional - use minBOF as objective
minBOFind = find(strcmp('BIOMASSMinimal',PanGEM.rxns));
PanGEM.c = zeros(numel(PanGEM.rxns),1);
PanGEM.c(minBOFind) = 1;

[essentialRxns, results] = solveCS(PanGEM);
% also don't mess with photosynthesis and oxphos... so let's add those to
% the list
PSInd = find(strcmp('Photosynthesis',PanGEM.subSystems));
OxPhosInd = find(strcmp('Oxidative phosphorylation',PanGEM.subSystems));
essentialRxns2 = unique([essentialRxns, PSInd',OxPhosInd'])
essentialRxns = essentialRxns2;
%% Get reactions to remove for each strain

    
for i = 1:nStrains;
    rxnsToRemove{i} = getRxnsToRemove(PanGEM,Dat, strain{i});
end

%%%%%%%%%%%%%% START Optional %%%%%%%%%%%%%%%%%%
% % Some sanity checks on the reactions to be removed
% AllToBeRemoved = unique(vertcat(rxnsToRemove{:}));
% for i = 1:numel(rxnsToRemove)
%     if i == 1
%         conservedSet = rxnsToRemove{i}; 
%     end    
%     newSet = rxnsToRemove{i};
%     conservedSet = intersect(conservedSet,newSet);
%     % keep track of the size of conserved sets
%     nConserved(i) = numel(conservedSet);
% end
% %
% % conservedSet should be the genes present in the pangenome that are not
% % found in the set of 'high quality' genomes. In other words, the genes
% % only found in SAGs and Metagenome derived genomes. Let's check that that
% % is true:
% 
% LowQualStrains = Dat.orgDatabase.Organism_Name(Dat.orgDatabase.Scaffolds>=30);
% for i = 1:numel(conservedSet)
% clear conservedKO conservedKO2 conservedKO3 conservedKOSplit conservedKOStrain isSplit isSplit2 isSplit3 isSAG
% 
%     conservedRxnInd = find(strcmp(conservedSet{i},PanGEM.rxns));
%     conservedKO = PanGEM.grRules(conservedRxnInd);
%     isSplit = strfind(conservedKO,'and');
%     isSplit2 = strfind(conservedKO,'or');
%     isSplit3 = ~isempty(isSplit{1}) | ~isempty(isSplit2{1});
%     if isSplit3>0
%         conservedKO2 = strrep(conservedKO,'(','');
%         conservedKO3 = strrep(conservedKO2,')','');
%         if ~isempty(isSplit{1})
%             conservedKOSplit = strsplit(conservedKO3{1},{' and '});
%         end
%         if ~isempty(isSplit2{1})
%             conservedKOSplit = strsplit(conservedKO3{1},{' or '});
%         end
%         
%         nConservedKO = numel(conservedKOSplit);
%         conservedKO = conservedKOSplit;
%     else nConservedKO = 1;
%     end
%     for j = 1:nConservedKO
%         tempInd = find(strcmp(conservedKO{j},Dat.KO));
%         conservedKOStrain = Dat.Strain(tempInd(1));
%         isSAG(j) = sum(strcmp(conservedKOStrain,LowQualStrains));
%     end
%     sumIsSAG{i} = sum(isSAG);    
% end
%   
% nRxnsToRemove = cellfun(@numel, rxnsToRemove);
% nGenes = Dat.orgDatabase.Genes(Dat.orgDatabase.Scaffolds<30);
% 
% figure
% h = plot(nGenes, nRxnsToRemove,'.k','MarkerSize',20);
% hold on
% h2 = text(nGenes,nRxnsToRemove,strrep(Dat.orgDatabase.Strain(Dat.orgDatabase.Scaffolds<30),'_',' '))
% xlabel('Number of Genes','FontSize',20);
% ylabel('Number of reactions to remove','FontSize',20);
% set(gca,'FontSize',20);
%%%%%%%%%%%%%% END Optional %%%%%%%%%%%%%%%%%%


%% Exclude essential reactions from reactions to remove and generate strain models
for i = 1:nStrains
    tempRxnstoRemove = rxnsToRemove{i};
    nRTR(i) = numel(tempRxnstoRemove);
    canRemove{i} = setdiff(tempRxnstoRemove,PanGEM.rxns(essentialRxns));
    cannotRemove{i} = intersect(tempRxnstoRemove,PanGEM.rxns(essentialRxns));
    pctRemoved(i) = 100*(numel(canRemove{i})/nRTR(i));
end


% let's have a look at the conserved set first
allCannotRemove = unique(vertcat(cannotRemove{:}))
coreCannotRemove{1} = cannotRemove{1};
for i = 1:nStrains-1
    coreCannotRemove{i+1} = intersect(coreCannotRemove{i},cannotRemove{i+1});
end
coreCannotRemove2 = coreCannotRemove{end};
% short list (13 total), so then each model has a different set of
% reactions which cannot be removed. Hmm...



% Populate strain models
% we'll do this two ways...
% 1. create reduced models (removing reactions)
% 2. by editing S in the PanGEM;

for i = 1:nStrains
    targStr = strain{i};
    redMod = removeReactions(PanGEM,canRemove{i},true,true,false);
    redMod.id = targStr;
    targStrMod2.(targStr) = redMod;
end

for i = 1:nStrains
    nRxnsFinal2(i) = numel(targStrMod2.(strain{i}).rxns);
end


for i = 1:nStrains
    targStr = strain{i};
    [junk, removeInd] = intersect(PanGEM.rxns,canRemove{i});
    targStrMod3.(targStr) = PanGEM;
    targStrMod3.(targStr).S(:,removeInd) = 0;
end

save('CBIOMES/Pangenomes/Prochlorococcus/StrMods_v20200105.mat','targStrMod3');
%% let's try paring back this number further using a random search

nPermutations = 50;
minGrowth = 0.1;
targStrMod4 = struct;
for i = 1:numel(strain)
    clear rxnsToRemoveResidual
    rxnsToRemoveResidual = setdiff(cannotRemove{i},coreCannotRemove2);
    [targStrMod4.(strain{i}),nRxns{i}] = getStrainMod_RandomSearch(targStrMod3.(strain{i}),rxnsToRemoveResidual,nPermutations,minGrowth);
end    

% check progress
nRxns2 = reshape(cell2mat(nRxns),50,46);
nRxns2(nRxns2==100000) = NaN;

imagesc(nRxns2)
colorbar

save('CBIOMES/Pangenomes/Prochlorococcus/StrMods_plus_v20190922.mat','targStrMod4');

    
%% Growth potential of each strain
nStrains = numel(strain);
tempMod = PanGEM; % add the whole PanGEM to the mix
[magicGrowth{1}, Carbon{1}, Nitrogen{1}, Phosphorus{1}, Sulfur{1}] = growthPotential(tempMod);

for i =2:nStrains +1;
    tempMod = targStrMod3.(strain{i-1});
    [magicGrowth{i}, Carbon{i}, Nitrogen{i}, Phosphorus{i}, Sulfur{i}] = growthPotential(tempMod);
end
nStrains = nStrains + 1; 
% kind of a pain to organize... here we go

for i = 1:nStrains
    temp{i} = Carbon{i}.substrates;
end
allCarbon = unique([vertcat(temp{:})]);

Carbon2 = zeros(numel(allCarbon),nStrains);
for i = 1:numel(allCarbon)
    for j = 1:nStrains
        tempInd = find(strcmp(allCarbon{i},Carbon{j}.substrates));
        if isempty(tempInd)
            Carbon2(i,j) = NaN;
        else
            Carbon2(i,j) = Carbon{j}.growth(tempInd);
        end
    end
end


for i = 1:nStrains
    temp{i} = Nitrogen{i}.substrates;
end
allNitrogen = unique([vertcat(temp{:})]);

Nitrogen2 = zeros(numel(allNitrogen),nStrains);
for i = 1:numel(allNitrogen)
    for j = 1:nStrains
        tempInd = find(strcmp(allNitrogen{i},Nitrogen{j}.substrates));
        if isempty(tempInd)
            Nitrogen2(i,j) = NaN;
        else
            Nitrogen2(i,j) = Nitrogen{j}.growth(tempInd);
        end
    end
end


for i = 1:nStrains
    temp{i} = Phosphorus{i}.substrates;
end
allPhosphorus = unique([vertcat(temp{:})]);

Phosphorus2 = zeros(numel(allPhosphorus),nStrains);
for i = 1:numel(allPhosphorus)
    for j = 1:nStrains
        tempInd = find(strcmp(allPhosphorus{i},Phosphorus{j}.substrates));
        if isempty(tempInd)
            Phosphorus2(i,j) = NaN;
        else
            Phosphorus2(i,j) = Phosphorus{j}.growth(tempInd);
        end
    end
end


for i = 1:nStrains
    temp{i} = Sulfur{i}.substrates;
end
allSulfur = unique([vertcat(temp{:})]);

Sulfur2 = zeros(numel(allSulfur),nStrains);
for i = 1:numel(allSulfur)
    for j = 1:nStrains
        tempInd = find(strcmp(allSulfur{i},Sulfur{j}.substrates));
        if isempty(tempInd)
            Sulfur2(i,j) = NaN;
        else
            Sulfur2(i,j) = Sulfur{j}.growth(tempInd);
        end
    end
end

AllGrowth = [Carbon2;Nitrogen2;Phosphorus2;Sulfur2]; % gDW biomass mol-1 h-1
AllSubstrates = [allCarbon;allNitrogen;allPhosphorus;allSulfur];
strainNames2= strrep(strain,'_','-');
strainNames3 = [{'PanGEM'}; strainNames2];
strainNames3 = strrep(strainNames3,'CCMP1375-SS120','SS120');
% Plot results

AllGrowthIO = zeros(size(AllGrowth));
AllGrowthIO(find(abs(AllGrowth)>.05)) = 1;
% let's just plot those substrates which support growth in the PanGEM
posGrowthInd = find(AllGrowthIO(:,1));
% let's also re-order into ecotypes
ecotypes = [{'LLII/III'},{'LLII/III'},{'LLII/III'},{'HLI'},{'LLIV'},...
    {'HLII'},{'LLI'},{'LLI'},{'HLII'},{'HLII'},{'HLII'},{'HLII'},{'HLI'},...
    {'HLII'},{'LLII/III'},{'HLII'},{'HLII'},{'HLII'},{'LLI'},{'LLII/III'},...
    {'HLII'},{'LLII/III'},{'LLIV'},{'HLII'},{'HLII'},{'HLI'},{'HLII'},...
    {'HLII'},{'HLII'},{'HLII'},{'HLII'},{'HLII'},{'HLII'},{'HLII'},...
    {'HLII'},{'LLII/III'},{'LLI'},{'LLII/III'},{'LLII/III'},{'LLII/III'},...
    {'LLIV'},{'LLIV'},{'LLIV'},{'LLIV'},{'LLIV'},{'LLIV'}];
ecoTemp = [{'HLI'},{'HLII'},{'LLI'},{'LLII/III'},{'LLIV'}];
for i = 1:numel(ecoTemp)
    group{i} = find(strcmp(ecoTemp{i},ecotypes))
end
group2 = vertcat([group{:}])
group3 = group2 + 1;
group4 = [1 group3];

Substrates2 = [{'Bicarb'},{'L-Alanine'},{'Ammonium'},{'L-Asparagine'},{'L-Aspartate'},{'Cyanate'},{'L-Glycine'},{'Nitrate'},{'Nitrite'},{'L-Serine'},{'Urea'},{'Urea2'},{'2-AEP'},{'Orthophosphate'},{'Phospholipids'},{'Phosphite'},{'Sulfate'},{'Sulfate2'},{'Thiosulfate'}]
SubstrateOrder = [3,8,9,11,6,2,4,5,7,10,14,16,15,17,19];
figure
imagesc(AllGrowthIO(posGrowthInd(SubstrateOrder),group4));
%title('Growth/no Growth')
set(gca,'YTick',1:numel(Substrates2(SubstrateOrder)),'YTickLabels',Substrates2(SubstrateOrder))
set(gca,'XTick',1:nStrains,'XTickLabels',strainNames3(group4),'XTickLabelRotation',90)
set(gca,'FontSize',20)
grid on
colormap('gray')




%% Analyze strain models BOF capabilities
BOFInd = find(strcmp('BIOMASSCOMBINED',PanGEM.rxns))
biometsInd = find(PanGEM.S(:,BOFInd));
biomets = PanGEM.mets(biometsInd);

for i = 1:numel(strain)
    tempMod = targStrMod3.(strain{i});
    
    % get model sizes
    modSize{i} = size(tempMod.S);
    
    % see if model can produce flexBOF mets
    for j = 1:numel(biomets)
        isPresent = find(strcmp(biomets{j},tempMod.mets));
        if ~isempty(isPresent)
            tempRes(j) = canProduce(tempMod,biomets{j});
            if tempRes(j)
                PanGEM_MetInd = find(strcmp(biomets{j},PanGEM.mets));
                coef = PanGEM.S(PanGEM_MetInd,BOFInd);
                tempMod_minBOFInd = find(strcmp('BIOMASSMinimal',tempMod.rxns));
                tempMod_flexBOF_MetInd = find(strcmp(biomets{j},tempMod.mets));
                tempMod.S(tempMod_flexBOF_MetInd,tempMod_minBOFInd) = coef;
            end
        end
    end
    % get a list of flexBOF metabolites that can be produced
    modFlexBOFMets{i} = biomets(tempRes);
    
    % get GPP, NPP, PSII/PSI, PQ from fluxes
    tempSol = solveLP(tempMod)
    if tempSol.stat == 1
        CO2EX_Ind = find(strcmp('CO2EX',tempMod.rxns));
        HCO3EX_Ind = find(strcmp('HCO3EX',tempMod.rxns));
        RuBisCO_Ind = find(strcmp('R00024',tempMod.rxns));
        PSII_Ind = find(strcmp('R09503',tempMod.rxns));
        PSIIabs_Ind = find(strcmp('PSIIabs',tempMod.rxns));
        PSIabs_Ind = find(strcmp('PSIabs',tempMod.rxns));
        growth(i) = tempSol.f;
        GPP(i) = (tempSol.x(RuBisCO_Ind));
        NPP(i) = GPP(i) - (tempSol.x(CO2EX_Ind) + tempSol.x(HCO3EX_Ind));
        PQ(i) = tempSol.x(PSII_Ind) / GPP(i);
        PSII_PSI_ratio(i) = tempSol.x(PSIIabs_Ind) / tempSol.x(PSIabs_Ind);
    else
        GPP(i) = NaN; NPP(i) = NaN; PQ(i) = NaN; PSII_PSI_ratio(i) = NaN; growth(i) = NaN;
    end
    clear tempRes 
end

% plot a canProduce matrix for each model
flexBOF_MetID_all = unique(vertcat(modFlexBOFMets{:}));

canProduceMat = zeros(numel(strain),numel(flexBOF_MetID_all));

for i = 1:numel(strain)
    for j = 1:numel(modFlexBOFMets{i})
        tempInd = find(strcmp(modFlexBOFMets{i}(j),flexBOF_MetID_all));
        if ~isempty(tempInd)
            canProduceMat(i,tempInd) = 1;
        else
            canProduceMat(i,tempInd) = 0;
        end
    end
end
canProduceMat = canProduceMat';
% add a vector of ones for PanGEM
canProduceMat2 = [ones(size(canProduceMat,1),1), canProduceMat];

figure
h = imagesc(canProduceMat2)
set(gca,'YTick',1:size(canProduceMat2,1),'YTickLabels',flexBOF_MetID_all)
set(gca,'XTick',1:nStrains,'XTickLabels',strainNames3(group4),'XTickLabelRotation',90)
set(gca,'FontSize',10)
grid on
colormap('gray')



% plot just the B vitamins
Vitamins = [{'B_1 Thiamine'},{'B_6 Pyridoxal'},{'B_7 Biotin'},{'B_1_2 Cobalamin'}];
VitOrder = [88, 80, 14, 23];

VitDatOnly = canProduceMat2(VitOrder,group4)

figure
h = imagesc(VitDatOnly)
set(gca,'YTick',1:numel(Vitamins),'YTickLabels',Vitamins)
set(gca,'XTick',1:nStrains,'XTickLabels',strainNames3(group4),'XTickLabelRotation',90)
set(gca,'FontSize',20)
grid on
colormap('gray')

%% Analyze strain models growth using an edited BOF specific to each strain
% use canProduceMat to generate a new BOF
strMods = targStrMod3; % targStrMod targStrMod2, targStrMod3
for i = 1:nStrains-1
    clear canProdMets
    tempMod = strMods.(strain{i});
    tempModMinBOFInd = find(strcmp('BIOMASSMinimal',tempMod.rxns));
    tempMod.c = zeros(numel(tempMod.rxns),1);
    tempMod.c(tempModMinBOFInd) = 1;
    % get indices of metabolites which can be produced for this strain
    canProdMets = flexBOF_MetID_all(find(canProduceMat(:,i)));
    % get the metInd for each
    [junk, canProdInd] = intersect(tempMod.mets,canProdMets);
    % get the corresponding coefficients in the full BOF
    tempModFullBOFInd = find(strcmp('BIOMASSCOMBINED',tempMod.rxns));
    tempCoefs = full(tempMod.S(canProdInd,tempModFullBOFInd));
    tempMod.S(canProdInd,tempModMinBOFInd) = tempCoefs;
    tempSol = solveLP(tempMod,1); % min sum of fluxes
    if tempSol.stat ==1;
        growth(i) = tempSol.f;
        fluxes{i} = tempSol.x;
    else growth(i) = 0;
        fluxes{i} = zeros(numel(tempMod.rxns),1);
    end
end



fluxes2 = cell2mat(fluxes);
fullSol = solveLP(PanGEM,1);
fluxDev = fluxes2 - repmat(fullSol.x,1,nStrains-1);
minflux = repmat(min(fluxDev')',1,size(fluxDev,2));
maxflux = repmat(max(fluxDev')',1,size(fluxDev,2));
scoreFlux = (fluxDev - minflux) ./ (maxflux - minflux);
nzFluxes = find(sum(scoreFlux,2)>1e-6)
figure
imagesc(scoreFlux(datasample(nzFluxes,numel(nzFluxes),'Replace',false),group2));
colorbar
grid on
set(gca,'XTick',1:nStrains-1,'XTickLabels',strainNames2(group2),'XTickLabelRotation',90)

% some specific fluxes and flux ratios
LightInd = find(strcmp('LightTRANS',PanGEM.rxns));
PSIabsInd = find(strcmp('PSIabs',PanGEM.rxns));
PSIIabsInd = find(strcmp('PSIIabs',PanGEM.rxns));
CO2TransInd = find(strcmp('CO2TRANS',PanGEM.rxns));
HCO3TransInd = find(strcmp('HCO3TRANScyt',PanGEM.rxns));
O2TransInd = find(strcmp('O2TRANS',PanGEM.rxns));
tempInd = find(strcmp('OrthophosphateTRANS',PanGEM.rxns));
figure
subplot(4,1,1)
h1 = bar(fluxes2(tempInd,:)./growth)

subplot(4,1,2)
h2 = bar(fluxes2(PSIabsInd,:) ./ fluxes2(PSIIabsInd,:))

subplot(4,1,3)
h3 = bar(fluxes2(O2TransInd,:) ./ ( fluxes2(CO2TransInd,:) + fluxes2(HCO3TransInd,:) ))

subplot(4,1,4)
h4 = bar(growth)



%% Analyze strain models growth using an edited BOF specific to each strain with a limitation
% use canProduceMat to generate a new BOF
strMods = targStrMod3; % targStrMod targStrMod2, targStrMod3
% Let's limit ammonia uptake
NH3TransInd = find(strcmp('AmmoniaTRANS',PanGEM.rxns));
NH3Uptake = sol.x(NH3TransInd);


for i = 1:nStrains-1
    clear canProdMets
    tempMod = strMods.(strain{i});
    tempModMinBOFInd = find(strcmp('BIOMASSMinimal',tempMod.rxns));
    tempMod.c = zeros(numel(tempMod.rxns),1);
    tempMod.c(tempModMinBOFInd) = 1;
    % get indices of metabolites which can be produced for this strain
    canProdMets = flexBOF_MetID_all(find(canProduceMat(:,i)));
    % get the metInd for each
    [junk, canProdInd] = intersect(tempMod.mets,canProdMets);
    % get the corresponding coefficients in the full BOF
    tempModFullBOFInd = find(strcmp('BIOMASSCOMBINED',tempMod.rxns));
    tempCoefs = full(tempMod.S(canProdInd,tempModFullBOFInd));
    tempMod.S(canProdInd,tempModMinBOFInd) = tempCoefs;
    % set limitation
    tempMod.lb(NH3TransInd) = NH3Uptake*.5;
    tempSol = solveLP(tempMod,1); % min sum of fluxes
    if tempSol.stat ==1;
        growth(i) = tempSol.f;
        fluxes{i} = tempSol.x;
    else growth(i) = 0;
        fluxes{i} = zeros(numel(tempMod.rxns),1);
    end
end



fluxes2 = cell2mat(fluxes);
fullSol = solveLP(PanGEM,1);
fluxDev = fluxes2 - repmat(fullSol.x,1,nStrains-1);
minflux = repmat(min(fluxDev')',1,size(fluxDev,2));
maxflux = repmat(max(fluxDev')',1,size(fluxDev,2));
scoreFlux = (fluxDev - minflux) ./ (maxflux - minflux);
nzFluxes = find(sum(scoreFlux,2)>1e-6)
figure
imagesc(scoreFlux(datasample(nzFluxes,numel(nzFluxes),'Replace',false),group2));
colorbar
grid on
set(gca,'XTick',1:nStrains-1,'XTickLabels',strainNames2(group2),'XTickLabelRotation',90)

% some specific fluxes and flux ratios
LightInd = find(strcmp('LightTRANS',PanGEM.rxns));
PSIabsInd = find(strcmp('PSIabs',PanGEM.rxns));
PSIIabsInd = find(strcmp('PSIIabs',PanGEM.rxns));
CO2TransInd = find(strcmp('CO2TRANS',PanGEM.rxns));
HCO3TransInd = find(strcmp('HCO3TRANScyt',PanGEM.rxns));
O2TransInd = find(strcmp('O2TRANS',PanGEM.rxns));
tempInd = find(strcmp('OrthophosphateTRANS',PanGEM.rxns));
figure
subplot(4,1,1)
h1 = bar(fluxes2(tempInd,:)./growth)

subplot(4,1,2)
h2 = bar(fluxes2(PSIabsInd,:) ./ fluxes2(PSIIabsInd,:))

subplot(4,1,3)
h3 = bar(fluxes2(O2TransInd,:) ./ ( fluxes2(CO2TransInd,:) + fluxes2(HCO3TransInd,:) ))

subplot(4,1,4)
h4 = bar(growth)





%%
% new order from getPhylogeneticTree.m
for i = 1:numel(leafNames)
    temp = find(strcmp(leafNames{i},strain));
    if ~isempty(temp)
        strainInd(i) = temp;
    else strainInd(i) = NaN;
    end
end
justPE = flexBOF_MetID_all(2:3);
canProduceMat2 = zeros(numel(strainInd),numel(justPE));
for i = 1:numel(strainInd);
    if ~isnan(strainInd(i))
        canProduceMat2(i,:) = canProduceMat(strainInd(i),2:3);
    else canProduceMat2(i,:) = repmat(-1,1,numel(justPE))
    end
end




h = plot(UPGMAtree);
pos = get(gcf, 'Position'); %// gives x left, y bottom, width, height
width = pos(3);
height = pos(4);
figure('Position',pos)
% subplot(1,5,1)
% genomeSize = Dat.orgDatabase.Size_Mb(strainInd);
% barh(genomeSize);
% set(gca,'YTick',1:numel(leafNames),'YTickLabel',strrep(leafNames,'_',' '))
% set(gca,'YDir','reverse')

% subplot(1,5,2)
% h = imagesc(canProduceMat2);
% xlab = strrep(justPE,'_',' ');
% ylab = strrep(leafNames,'_',' ');
% set(gca,'XTick',1:numel(xlab),'XTickLabels',xlab)
% set(gca,'YTick',1:numel(ylab),'YTickLabels',ylab)
% xtickangle(90)

subplot(1,3,1)
barh(nRxns(strainInd));
%set(gca,'YTick',1:numel(leafNames),'YTickLabel',strrep(leafNames,'_',' '))
set(gca,'XLim',[850 1100]);
set(gca,'YDir','reverse')
set(gca,'YTick',[])


subplot(1,3,2)

barh(nMets(strainInd));
%set(gca,'YTick',1:numel(leafNames),'YTickLabel',strrep(leafNames,'_',' '))
set(gca,'XLim',[800 1000]);
set(gca,'YDir','reverse')
set(gca,'YTick',[])

subplot(1,3,3)
barh(NPP(strainInd));
set(gca,'YTick',1:numel(leafNames),'YTickLabel',strrep(leafNames,'_',' '))
set(gca,'XLim',[8.57 8.64]);
set(gca,'YDir','reverse')
%set(gca,'XTick',[],'YTick',[])

tightfig



%% Elementary flux modes
modStr = struct;
modStr.stoich = full(PanGEM.S);
modStr.reversibilities = PanGEM.ub == 0 | PanGEM.lb == 0;
modStr.metaboliteNames = PanGEM.mets;
modStr.reactionNames = PanGEM.rxns;
rmpath Software/RAVEN-2.0.0-alpha/software/apache-poi/
cd Software/efmtool
mnet = CalculateFluxModes(modStr)
cd ../../
addpath Software/RAVEN-2.0.0-alpha/software/apache-poi/


% ran out of memory

