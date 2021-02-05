%% Manual curation part 2
% finding blocked reactions, canProduce/consume, finding gaps, filling, etc


%% solve

fileName = sprintf('%s','CBIOMES/Pangenomes/Prochlorococcus/PanGEM_pro_20200105.xlsx');
PanGEM = importExcelModel(fileName,true,true,false);
%tic
%[sol hsSolOut] = solveLP(PanGEM);
%t1 = toc
[sol hsSolOut] = solveLP(PanGEM,1)
%t2 = toc

condMod = PanGEM;
condMod.ub(find(strcmp('R00475',condMod.rxns))) = 0;
condMod.lb(find(strcmp('LightEX',condMod.rxns))) = -10;

condSol = solveLP(condMod,1);

%% map Reactions
map={'map01040'};
for i=1:numel(map)
    drawKEGGPathwayKO(PanGEM,cell2mat(map(i)));
end

map={'map00770'};
for i=1:numel(map)
    drawKEGGPathwayKO(targStrMod3.(strain{7}),cell2mat(map(i)));
end
%% map RN's
map={'map00270'};
for i=1:numel(map)
    drawKEGGPathwayKO(PanGEM,cell2mat(map(i)));
end
%% canProduce/consume
clear notProduced notConsumed
produced = canProduce(PanGEM); 
consumed = canConsume(PanGEM);
notProduced = PanGEM.mets(~produced);
notConsumed = PanGEM.mets(~consumed);
notProducedConsumed = intersect(notProduced,notConsumed)
% lets see how many of the biomets we can produce
BOFInd = find(strcmp('BIOMASSCRUDE',PanGEM.rxns))
biometsInd = find(PanGEM.S(:,BOFInd));
biomets = PanGEM.mets(biometsInd);
biometsNotProduced = intersect(biomets,notProduced)
biometsNotConsumed = intersect(biomets,notConsumed)

% use haveFlux to detect reactions which cannot carry flux
J=haveFlux(PanGEM,10^-40);
Flux = [PanGEM.rxns(J) PanGEM.rxnNames(J) constructEquations(PanGEM,PanGEM.rxns(J)) PanGEM.eccodes(J) PanGEM.subSystems(J)];
noFlux = [PanGEM.rxns(~J) PanGEM.rxnNames(~J) constructEquations(PanGEM,PanGEM.rxns(~J)) PanGEM.eccodes(~J) PanGEM.subSystems(~J)];


% let's see if any of the biomets are holding up growth
for i = 1:numel(biomets);
    tempMod = PanGEM;
    tempMod.S(biometsInd(i),BOFInd) = 0;
    tempSol = solveLP(tempMod);
    tempGrowth(i) = tempSol.f;
end

% cool, one hit! let's see what it is:
PanGEM.mets(biometsInd(find(tempGrowth<-1)))
% thiamin diphosphate again! sonofa
% look for disconnected reactins
for i = 1:numel(PanGEM.rxns);
    rxnName = PanGEM.rxns{i};
    tempStruct = checkRxn(PanGEM,rxnName,10e-40);
    nReactants = numel(tempStruct.reactants);
    nProducts = numel(tempStruct.products);
    nCanMake = sum(tempStruct.canMake);
    nCanConsume = sum(tempStruct.canConsume);
    if nReactants == nCanMake & nProducts == nCanConsume;
        goodRxn(i) = 1;
    else goodRxn(i) = 0;
    end
end
goodRxn = goodRxn';


% run a gap report
[noFluxRxns, noFluxRxnsRelaxed, subGraphs, notProducedMets, minToConnect,...
    neededForProductionMat, canProduceWithoutInput, canConsumeWithoutOutput, ...
    connectedFromTemplates, addedFromTemplates]=gapReport(PanGEM);
% takes a while so lets save this...
gapReport_PanGEM = struct;
gapReport_PanGEM.noFluxRxns = noFluxRxns;
gapReport_PanGEM.noFluxRxnsRelaxed = noFluxRxnsRelaxed;
gapReport_PanGEM.subGraphs = subGraphs;
gapReport_PanGEM.notProducedMets = notProducedMets;
gapReport_PanGEM.minToConnect = minToConnect;
gapReport_PanGEM.neededForProductionMat = neededForProductionMat;
save('CBIOMES/Pangenomes/SAR86/gapReport.mat','gapReport_PanGEM');





%%
% now let's set priorities, we can do this by iteratively adding exchange
% reactions for metabolites that cannot be produced, and then checking for
% the change in the number of metabolites that can subsequently be
% produced. We'll rank this list in order to focus on where the big
% problems are.


for i = 1:numel(notProduced);
    tempMod = PanGEM;
    tempMod = addExchangeRxns(tempMod,'in',notProduced{i});
    tempProduced = canProduce(tempMod);
    nProduced(i) = sum(tempProduced);
end



%% Extras
load('Pro_Assembly3.mat');


ind = find(strcmp(Dat.KO,'K15577'));[Dat.GeneID(ind) Dat.GeneVerbose(ind) Dat.Strain(ind)]
%% Stop the leaks!
% Find where illegal carbon fixation is going on
CO2Ind = find(strcmp('CO2',PanGEM.mets));
CO2RxnsInd = find(PanGEM.S(CO2Ind,:));
CO2RxnsCoef = full(PanGEM.S(CO2Ind,:));
CO2RxnsCoef2 = CO2RxnsCoef(CO2RxnsInd);

CO2ConsumptionRxnsInd = intersect(find(CO2RxnsCoef2<0),find(sol.x(CO2RxnsInd)>0))
CO2ConsumptionRxnsInd2 = intersect(find(CO2RxnsCoef2>0),find(sol.x(CO2RxnsInd)<0))
CO2ConsumptionRxnsInd3 = [CO2ConsumptionRxnsInd; CO2ConsumptionRxnsInd2]

a = PanGEM.rxns(CO2RxnsInd(CO2ConsumptionRxnsInd3));
b = sol.x(CO2RxnsInd(CO2ConsumptionRxnsInd3));

table(a,b)

%% Change rxn bounds to find out which reactions are suspect

for i = 1:numel(PanGEM.rxns);
    tempMod = PanGEM;
    if tempMod.lb(i) == 0;
        tempMod.lb(i) = -1000;
        tempSol = solveLP(tempMod);
        if tempSol.stat == 1;
            lbSol(i) = tempSol.f;
        else lbSol(i) = 0;
        end
    else lbSol(i) = sol.f;
    end
end

for i = 1:numel(PanGEM.rxns);
    tempMod = PanGEM;
    if tempMod.ub(i) == 0;
        tempMod.ub(i) = 1000;
        tempSol = solveLP(tempMod);
        if tempSol.stat == 1;
            ubSol(i) = tempSol.f;
        else ubSol(i) = 0;
        end
    else ubSol(i) = sol.f;
    end
end


[crap1, LBSuspect] = sort(abs(lbSol),'descend')
[crap2, UBSuspect] = sort(abs(ubSol),'descend')

table(PanGEM.rxns(LBSuspect(1:50)),crap1(1:50)')
table(PanGEM.rxns(UBSuspect(1:50)),crap2(1:50)')
%% Check shadow prices to see what is limiting

[crap shadowInd] = sort(abs(hsSolOut.y),'descend')

table(PanGEM.mets(shadowInd(1:50)),crap(1:50))

%% Check where ATP is being produced
ATPInd = find(strcmp('ATP',PanGEM.mets));
ATPRxnsInd = find(PanGEM.S(ATPInd,:));
ATPRxnsCoef = full(PanGEM.S(ATPInd,:));
ATPRxnsCoef2 = ATPRxnsCoef(ATPRxnsInd);

ATPProductionRxnsInd = intersect(find(ATPRxnsCoef2>0),find(sol.x(ATPRxnsInd)>0))
ATPProductionRxnsInd2 = intersect(find(ATPRxnsCoef2<0),find(sol.x(ATPRxnsInd)<0))
ATPProductionRxnsInd3 = [ATPProductionRxnsInd; ATPProductionRxnsInd2]

a = PanGEM.rxns(ATPRxnsInd(ATPProductionRxnsInd3));
b = sol.x(ATPRxnsInd(ATPProductionRxnsInd3));

[bSort bInd] = sort(b,'descend');
aSort = a(bInd);

table(aSort,bSort)


%% Let's compare iJC568 and the PanGEM rxn directionality


fileName = sprintf('%s','Chalmers/Datasets/MED4_models/draftMED4/editing/draftMED4_20161114Combined.xlsx')
%fileName = sprintf('%s','Chalmers/Datasets/MED4_models/draftMED4/editing/iJC581_Excel.xlsx')
iJC568 = importExcelModel(fileName,true,true,false);
unused = iJC568.mets(~any(iJC568.S,2));
iJC568 = removeMets(iJC568,unused);
clear matchVec
for i = 1:numel(PanGEM.rxns)
    matchInd = find(strcmp(PanGEM.rxns{i},iJC568.rxns));
    if isempty(matchInd);
        matchVec{i} = 'Not Present';
    elseif PanGEM.lb(i) == iJC568.lb(matchInd) && PanGEM.ub(i) == iJC568.ub(matchInd)
        matchVec{i} = 'Match';
    else matchVec{i} = 'No Match';
    end
end

matchVec = matchVec'
    
    
%% YEAAAAAAAAAAAAH!!!

    

%% Check a list of reactions
% paste a list from excel into an empty cell array called temp
temp = {};

for i = 1:numel(temp)
    % check if can carry flux
    report = checkRxn(model,temp{i},10^-10);
    if sum(report.canMake) == numel(report.canMake);
        canMake(i) = 1;
    else canMake(i) = 0;
    end
    if sum(report.canConsume) == numel(report.canConsume);
        canConsume(i) = 1;
    else canConsume(i) = 0;
    end
end

    

%% clean up unused metabolites
% get a logical vector of whether a metabolite is to be removed from the
% model.
temp = ~any(full(PanGEM.S),2)



%% Get reaction directions from eQuillibrator database
model = PanGEM;

% exclude some reactions which were manually curated

% lets try with only those reactions which don't produce/consume cofactors
cofactors = [{'ATP'},{'ADP'},{'AMP'},{'NAD'},{'NADH'},{'NADP'},{'NADPH'},{'S_Adenosyl_L_methionine'},{'FAD'}];
for i = 1:numel(cofactors)
    cofactors_ind(i) = find(strcmp(cofactors{i},model.mets));
end
[junk, cofactors_rxnInd] = find(full(model.S(cofactors_ind,:)));
% lets also not mess with oxygen and co2
cofactors2 = [{'CO2'},{'Oxygen'}];
for i = 1:numel(cofactors2)
    cofactors_ind2(i) = find(strcmp(cofactors2{i},model.mets));
end
[junk, cofactors_rxnInd2] = find(full(model.S(cofactors_ind2,:)));
% and don't mess with photosynthesis and oxphos
PS_RxnInd = find(strcmp('Photosynthesis',model.subSystems));
OxPhos_RxnInd = find(strcmp('Oxidative phosphorylation',model.subSystems));

excluded = unique([cofactors_rxnInd; cofactors_rxnInd2; PS_RxnInd; OxPhos_RxnInd]);

cutoff = 10;
[lb, ub, change_bounds_ind] = getThermoDirectionality(model,excluded, cutoff);

model_new = model;
model_new.lb(change_bounds_ind) = -1000;
model_new.ub(change_bounds_ind) = 1000;

sol_new = solveLP(model_new);
sol_new2 = solveLP(model_new,1);

% work in excel
temp = zeros(numel(model_new.rxns),1);
temp(change_bounds_ind) = 1;




%% Check tasks
inputFile = 'CBIOMES/Scripts/ManualCuration/CheckTasksProPanGEM.xlsx';
taskStruct=parseTaskList(inputFile);

[taskReport essentialRxns taskStructure]=checkTasks(PanGEM,[],true,false,true,taskStruct)





%% Fix the light issue
% an issue with ATP and NAD(P)H production occurred when limiting light...
% let's hunt it down

fileName = sprintf('%s','CBIOMES/Pangenomes/Prochlorococcus/PanGEM_pro_20180205.xlsx');
PanGEM = importExcelModel(fileName,true,true,false);
sol = solveLP(PanGEM,1);

LightTransInd = find(strcmp('LightTRANS',PanGEM.rxns));
LightUptake = sol.x(LightTransInd);


LightSeries = 0:LightUptake/10:LightUptake + LightUptake/10;
tempMod = PanGEM;
for i = 1:numel(LightSeries);
    tempMod.ub(LightTransInd) = LightSeries(i);
    tempSol{i} = solveLP(tempMod,1);
    if tempSol{i}.stat==1
        growth(i) = tempSol{i}.f;
    else
        growth(i) = 0;
    end
end

% plot exchange fluxes
exInd = find(strcmp('Exchange',PanGEM.subSystems));
tempMat = zeros(numel(exInd),numel(LightSeries));
for i = 1:numel(LightSeries);
    tempMat(:,i) = tempSol{i}.x(exInd);
end

% get rid of zero columns
tempMatSum = sum(tempMat,2);
nzRows = find(abs(tempMatSum)>1e-9);


figure
subplot(1,2,1)
imagesc(tempMat(nzRows,:));
colorbar
set(gca,'YTick',1:numel(nzRows),'YTickLabels',PanGEM.rxns(exInd(nzRows)))

subplot(1,2,2)
plot(LightSeries,growth)

% found one... it was R09691 in carotenoid biosynthesis which was allowing
% ferredoxin to be reduced. 


%% compare 2 versions of the model (no change in S)


fileName = sprintf('%s','CBIOMES/Pangenomes/Prochlorococcus/PanGEM_pro_20180315.xlsx');
PanGEMA = importExcelModel(fileName,true,true,false);
solA = solveLP(PanGEMA,1);

fileName = sprintf('%s','CBIOMES/Pangenomes/Prochlorococcus/PanGEM_pro_20180320.xlsx');
PanGEMB = importExcelModel(fileName,true,true,false);
solB = solveLP(PanGEMB,1);

figure
plot(solA.x,solB.x,'.k','MarkerSize',20)

% get the changed fluxes
diffInd = find(abs(solA.x - solB.x) > 1); % larger than 1 mmol gDW-1 h-1

PanGEMA.rxns(diffInd)









