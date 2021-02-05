function [magicGrowth, Carbon, Nitrogen, Phosphorus, Sulfur] = growthPotential(model,excludedRxns)
%% Growth potential
% identify and quantify growth on sole C,N,P, and S substrates

%model = PanGEM;

%% get transporter reactions for each element
transRxnInd = find(strcmp('Transport',model.subSystems));

% get transport mets
[transMetInd, junk] = find(full(model.S(:,transRxnInd)));
% screen out the intracellular mets
extracellularInd = find(model.metComps == 3);
transMetInd2 = intersect(transMetInd,extracellularInd);

%% Get excluded reactions indices
excludedRxnsInd = find(contains(model.rxns,excludedRxns));
transMets = model.mets(transMetInd2);

% parse out formulas
formulas = model.metFormulas(transMetInd2);
[elements, useMat, exitFlag, MW] = parseFormulas(formulas, false,false,true);

% get indicess of CNPS
bioElements = [{'C'},{'N'},{'P'},{'S'}];
for i = 1:numel(bioElements);
    bioElementsInd(i) = find(strcmp(bioElements{i},elements.abbrevs));
end

% CNPS sources
CInd = find(useMat(:,bioElementsInd(1)));
CSources = transMets(CInd);
Catoms = useMat(CInd,bioElementsInd(1));

NInd = find(useMat(:,bioElementsInd(2)));
NSources = transMets(NInd);
Natoms = useMat(NInd,bioElementsInd(2));

PInd = find(useMat(:,bioElementsInd(3)));
PSources = transMets(PInd);
Patoms = useMat(PInd,bioElementsInd(3));

SInd = find(useMat(:,bioElementsInd(4)));
SSources = transMets(SInd);
Satoms = useMat(SInd,bioElementsInd(4));

% now link back to their respective transporter reactions
for i = 1:numel(CSources)
    temp = find(strcmp(CSources{i},model.mets));
    CSourceRxnInd = find(full(model.S(temp,:)));
    nTemp = numel(CSourceRxnInd);
    % make sure it's a transporter
    CSourceRxnInd2{i} = intersect(CSourceRxnInd,transRxnInd);
    CSourceRxns{i} = model.rxns(CSourceRxnInd2{i});
    if nTemp ==1
        CatomRxns{i} = Catoms(i);
    else
    CatomRxns{i} = repmat(Catoms(i),nTemp-1,1);
    end
end
CSourceRxns2a = vertcat(CSourceRxns{:});
CSourceRxnInd3a = vertcat(CSourceRxnInd2{:});
CatomRxn2a = vertcat(CatomRxns{:});

for i = 1:numel(NSources)
    temp = find(strcmp(NSources{i},model.mets));
    NSourceRxnInd = find(full(model.S(temp,:)));
    nTemp = numel(NSourceRxnInd);
    % make sure it's a transporter
    NSourceRxnInd2{i} = intersect(NSourceRxnInd,transRxnInd);
    NSourceRxns{i} = model.rxns(NSourceRxnInd2{i});
    if nTemp ==1
        NatomRxns{i} = Natoms(i);
    else
    NatomRxns{i} = repmat(Natoms(i),nTemp-1,1);
    end
end
NSourceRxns2a = vertcat(NSourceRxns{:});
NSourceRxnInd3a = vertcat(NSourceRxnInd2{:});
NatomRxn2a = vertcat(NatomRxns{:});

for i = 1:numel(PSources)
    temp = find(strcmp(PSources{i},model.mets));
    PSourceRxnInd = find(full(model.S(temp,:)));
    nTemp = numel(PSourceRxnInd);
    % make sure it's a transporter
    PSourceRxnInd2{i} = intersect(PSourceRxnInd,transRxnInd);
    PSourceRxns{i} = model.rxns(PSourceRxnInd2{i});
    if nTemp ==1
        PatomRxns{i} = Patoms(i);
    else
    PatomRxns{i} = repmat(Patoms(i),nTemp-1,1);
    end
end
PSourceRxns2a = vertcat(PSourceRxns{:});
PSourceRxnInd3a = vertcat(PSourceRxnInd2{:});
PatomRxn2a = vertcat(PatomRxns{:});

for i = 1:numel(SSources)
    temp = find(strcmp(SSources{i},model.mets));
    SSourceRxnInd = find(full(model.S(temp,:)));
    nTemp = numel(SSourceRxnInd);
    % make sure it's a transporter
    SSourceRxnInd2{i} = intersect(SSourceRxnInd,transRxnInd);
    SSourceRxns{i} = model.rxns(SSourceRxnInd2{i});
    if nTemp ==1
        SatomRxns{i} = Satoms(i);
    else
    SatomRxns{i} = repmat(Satoms(i),nTemp-1,1);
    end
end
SSourceRxns2a = vertcat(SSourceRxns{:});
SSourceRxnInd3a = vertcat(SSourceRxnInd2{:});
SatomRxn2a = vertcat(SatomRxns{:});

%% Exclude required transport reactions reactions
CSourceRxns2 = setdiff(CSourceRxns2a,excludedRxns);
[CSourceRxnInd3, ia] = setdiff(CSourceRxnInd3a,excludedRxnsInd);
CatomRxn2 = CatomRxn2a(ia);

NSourceRxns2 = setdiff(NSourceRxns2a,excludedRxns);
[NSourceRxnInd3, ia] = setdiff(NSourceRxnInd3a,excludedRxnsInd);
NatomRxn2 = NatomRxn2a(ia);

PSourceRxns2 = setdiff(PSourceRxns2a,excludedRxns);
[PSourceRxnInd3, ia] = setdiff(PSourceRxnInd3a,excludedRxnsInd);
PatomRxn2 = PatomRxn2a(ia);

SSourceRxns2 = setdiff(SSourceRxns2a,excludedRxns);
[SSourceRxnInd3, ia] = setdiff(SSourceRxnInd3a,excludedRxnsInd);
SatomRxn2 = SatomRxn2a(ia);

%% Make sure that the model cannot grow without any source of each element
sol2 = solveLP(model,1);

% Growth on no carbon source
tempMod = model;
tempMod.lb(CSourceRxnInd3) = 0;
tempMod.ub(CSourceRxnInd3) = 0;
ClimMod = tempMod;
noCsol = solveLP(ClimMod,1);

% Growth on no nitrogen source
tempMod = model;
tempMod.lb(NSourceRxnInd3) = 0;
tempMod.ub(NSourceRxnInd3) = 0;
NlimMod = tempMod;
noNsol = solveLP(NlimMod,1);

% Growth on no phosphorus source
tempMod = model;
tempMod.lb(PSourceRxnInd3) = 0;
tempMod.ub(PSourceRxnInd3) = 0;
PlimMod = tempMod;
noPsol = solveLP(PlimMod,1);

% Growth on no sulfur source
tempMod = model;
tempMod.lb(SSourceRxnInd3) = 0;
tempMod.ub(SSourceRxnInd3) = 0;
SlimMod = tempMod;
noSsol = solveLP(SlimMod,1);

% all ok?
threshold = 1e-8;
if abs(noCsol.f) < threshold && abs(noNsol.f) < threshold && abs(noPsol.f) < threshold && abs(noSsol.f) < threshold
    magicGrowth = 'No'
else magicGrowth = 'Yes';
end

%% Growth on each substrate as sole element source

% Carbon sources
for i = 1:numel(CSourceRxnInd3)
    clear tempMetsInd tempMetsComp extracellularMetInd extracellularMetInd2 extracellularMetInd3 eMetS tempSol
    tempMod = ClimMod;
    % get the sense of each reaction
    tempMetsInd = find(full(tempMod.S(:,CSourceRxnInd3(i))));
    tempMetsComp = tempMod.metComps(tempMetsInd);
    extracellularMetInd = find(tempMetsComp==3);
    extracellularMetInd2 = tempMetsInd(extracellularMetInd);
    % sometimes there are symporters, let's retrieve the CNPS substrates
    % only
    extracellularMetInd3 = intersect(extracellularMetInd2,transMetInd2);
    eMetS = full(tempMod.S(extracellularMetInd3,CSourceRxnInd3(i)));
    if eMetS < 0
        tempMod.lb(CSourceRxnInd3(i)) = 0;
        tempMod.ub(CSourceRxnInd3(i)) = 1000;
    elseif eMetS > 0
        tempMod.lb(CSourceRxnInd3(i)) = -1000;
        tempMod.ub(CSourceRxnInd3(i)) = 0;
    end
    %tempMod.lb(CSourceRxnInd3(i)) = -1000;
    %tempMod.ub(CSourceRxnInd3(i)) = 1000;
    tempSol = solveLP(tempMod,1);
    if tempSol.stat ==1
        Cgrowth(i) = tempSol.f / CatomRxn2(i);
    else Cgrowth(i) = NaN;
    end
end

% Nitrogen sources
for i = 1:numel(NSourceRxnInd3)
    clear tempMetsInd tempMetsComp extracellularMetInd extracellularMetInd2 extracellularMetInd3 eMetS tempSol
    tempMod = NlimMod;
    % get the sense of each reaction
    tempMetsInd = find(full(tempMod.S(:,NSourceRxnInd3(i))));
    tempMetsComp = tempMod.metComps(tempMetsInd);
    extracellularMetInd = find(tempMetsComp==3);
    extracellularMetInd2 = tempMetsInd(extracellularMetInd);
    % sometimes there are symporters, let's retrieve the CNPS substrates
    % only
    extracellularMetInd3 = intersect(extracellularMetInd2,transMetInd2);
    eMetS = full(tempMod.S(extracellularMetInd3,NSourceRxnInd3(i)));
    if eMetS < 0
        tempMod.lb(NSourceRxnInd3(i)) = 0;
        tempMod.ub(NSourceRxnInd3(i)) = 1000;
    elseif eMetS > 0
        tempMod.lb(NSourceRxnInd3(i)) = -1000;
        tempMod.ub(NSourceRxnInd3(i)) = 0;
    end
    tempSol = solveLP(tempMod,1);
    if tempSol.stat ==1
        Ngrowth(i) = tempSol.f / NatomRxn2(i);
    else Ngrowth(i) = NaN;
    end
end


% Phosphorus sources
for i = 1:numel(PSourceRxnInd3)
    clear tempMetsInd tempMetsComp extracellularMetInd extracellularMetInd2 extracellularMetInd3 eMetS tempSol
    tempMod = PlimMod;
    % get the sense of each reaction
    tempMetsInd = find(full(tempMod.S(:,PSourceRxnInd3(i))));
    tempMetsComp = tempMod.metComps(tempMetsInd);
    extracellularMetInd = find(tempMetsComp==3);
    extracellularMetInd2 = tempMetsInd(extracellularMetInd);
    % sometimes there are symporters, let's retrieve the CNPS substrates
    % only
    extracellularMetInd3 = intersect(extracellularMetInd2,transMetInd2);
    eMetS = full(tempMod.S(extracellularMetInd3,PSourceRxnInd3(i)));
    if eMetS < 0
        tempMod.lb(PSourceRxnInd3(i)) = 0;
        tempMod.ub(PSourceRxnInd3(i)) = 1000;
    elseif eMetS > 0
        tempMod.lb(PSourceRxnInd3(i)) = -1000;
        tempMod.ub(PSourceRxnInd3(i)) = 0;
    end
    tempSol = solveLP(tempMod,1);
    if tempSol.stat ==1
        Pgrowth(i) = tempSol.f / PatomRxn2(i);
    else Pgrowth(i) = NaN;
    end
end


% Sulfur sources
for i = 1:numel(SSourceRxnInd3)
    clear tempMetsInd tempMetsComp extracellularMetInd extracellularMetInd2 extracellularMetInd3 eMetS tempSol
    tempMod = SlimMod;
    % get the sense of each reaction
    tempMetsInd = find(full(tempMod.S(:,SSourceRxnInd3(i))));
    tempMetsComp = tempMod.metComps(tempMetsInd);
    extracellularMetInd = find(tempMetsComp==3);
    extracellularMetInd2 = tempMetsInd(extracellularMetInd);
    % sometimes there are symporters, let's retrieve the CNPS substrates
    % only
    extracellularMetInd3 = intersect(extracellularMetInd2,transMetInd2);
    eMetS = full(tempMod.S(extracellularMetInd3,SSourceRxnInd3(i)));
    if eMetS < 0
        tempMod.lb(SSourceRxnInd3(i)) = 0;
        tempMod.ub(SSourceRxnInd3(i)) = 1000;
    elseif eMetS > 0
        tempMod.lb(SSourceRxnInd3(i)) = -1000;
        tempMod.ub(SSourceRxnInd3(i)) = 0;
    end
    tempSol = solveLP(tempMod,1);
    if tempSol.stat ==1
        Sgrowth(i) = tempSol.f / SatomRxn2(i);
    else Sgrowth(i) = NaN;
    end
end


%% output
Carbon = struct;
Carbon.substrates = CSourceRxns2;
Carbon.growth = Cgrowth';
Carbon.atoms = CatomRxn2;

Nitrogen = struct;
Nitrogen.substrates = NSourceRxns2;
Nitrogen.growth = Ngrowth';
Nitrogen.atoms = NatomRxn2;

Phosphorus = struct;
Phosphorus.substrates = PSourceRxns2;
Phosphorus.growth = Pgrowth';
Phosphorus.atoms = PatomRxn2;

Sulfur = struct;
Sulfur.substrates = SSourceRxns2;
Sulfur.growth = Sgrowth';
Sulfur.atoms = SatomRxn2;
