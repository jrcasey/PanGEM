%% PanGEM_fromAssembly_Pro

FileNames.L1_Assembly = 'data/assemblies/Assembly_L1_20200815.mat';
FileNames.PanGEM = 'data/models/Pro_PanGEM.mat';
FileNames.keepers = 'data/assemblies/KO_toAdd_20200818.csv';
FileNames.KO_remove = 'data/assemblies/KO_toRemove_20200819.csv';
%% Load L1 Assembly

load(FileNames.L1_Assembly);
strains = Pro_Assembly.orgDatabase.StrainName;
nStr = numel(strains);

%% Load PanGEM
PanGEM = importExcelModel('/Users/jrcasey/Documents/MATLAB/CBIOMES/Pangenomes/Prochlorococcus/PanGEM_pro_20200527.xlsx',true,true,false)
%load(FileNames.PanGEM);

% check for growth
tempSol = solveLP(PanGEM,1);
nNZ_fluxes = numel(find(tempSol.x~=0));
pctNZ_fluxes = nNZ_fluxes ./ numel(tempSol.x);

%% Retrieve KEGG API KO<->RN relational table

% call kegg's REST API
response = urlread('http://rest.kegg.jp/link/reaction/ko');

% parse response
response2 = textscan(response,'%s%s', 'Delimiter','\t');
KO_RN.KO = strrep(response2{1},'ko:','');
KO_RN.RN = strrep(response2{2},'rn:','');

%% Get all unique reactions that match to KO's in the KEGG database
for a = 1:numel(Pro_Assembly.uniqueKO)
    hit_idx = find(strcmp(Pro_Assembly.uniqueKO{a},KO_RN.KO));
    nHits(a) = numel(hit_idx);
    RN_all{a} = KO_RN.RN(hit_idx);
end

% get unique reactions from all hits
uniqueRN = unique(vertcat(RN_all{:}));
n_uniqueRN = numel(uniqueRN);

% KO's with no reaction hits (to be inspected manually)
KO_noHits = Pro_Assembly.uniqueKO(nHits==0);

% let's whittle this list down by keeping all those KO's with no hits to
% KEGG reactions, but with hits to PanGEM reactions
KO_noHits_keep = intersect(KO_noHits,PanGEM.genes);

% Not much help (only 143 meet these criteria), so let's manually browse
% through the annotations
KO_stillNoHits = setdiff(KO_noHits,KO_noHits_keep);
listKO_noHits = strjoin(KO_stillNoHits,'+');
for a = 1:numel(KO_stillNoHits)
    try tempResponse = urlread(strcat('http://rest.kegg.jp/list/',KO_stillNoHits{a}));
    % parse response
    tempResponse2 = textscan(tempResponse,'%s%s', 'Delimiter','\t');
    col1(a) = tempResponse2{1};
    col2(a) = tempResponse2{2};
    end
end

% some have EC numbers, let's parse those out as well
for a = 1:numel(col2)
    tempStr = col2{a};
    if ~isempty(tempStr)
        startIdx = regexp(tempStr,'[EC:');
        if ~isempty(startIdx)
            hasEC(a) = 1;
            endIdx = regexp(tempStr,']');
            EC_vec = tempStr(startIdx+4:endIdx(end)-1);
            EC_vec2 = strsplit(EC_vec,' ');
            nEC = numel(EC_vec2);
            for b = 1:nEC
                col3{a}{b} = EC_vec2{b};
                
                try tempResponse = urlread(strcat('http://rest.kegg.jp/link/rn/',col3{a}{b}));
                    tempResponse2 = textscan(tempResponse,'%s%s', 'Delimiter','\t');
                    tempResponse3 = strrep(tempResponse2{2},'rn:','');
                    col4{a}{b} = strjoin(tempResponse3,';');
                end
            end
        end
    end
end

% gotta loop back through and get rid of duplicates
for a =1:numel(col4)
    nEC = numel(col4{a});
    if nEC >1
        EC_vec = strjoin([col4{a}],';');
        % now split up again
        EC_vec2 = strsplit(EC_vec,';');
        EC_vec3 = unique(EC_vec2);
        col4{a} = EC_vec3;
    end
end




% now let's append KO's which have a matching RN in PanGEM to grRules
temp_grRules = PanGEM.grRules;
temp_grRules2 = strrep(temp_grRules,'(','');
temp_grRules3 = strrep(temp_grRules2,')','');
temp_grRules4 = strrep(temp_grRules3,' or ',';');
temp_grRules5 = temp_grRules4;

for a =1:numel(col4)
    tempRN = col4{a};
    if ~isempty(tempRN)
        tempRN2 = strsplit(tempRN{1},';');
        nRN = numel(tempRN2);
        for b = 1:nRN
            PanGEM_RNidx = find(strcmp(tempRN2{b},PanGEM.rxns));
            if ~isempty(PanGEM_RNidx)
                PanGEM_KO_idx = find(PanGEM.rxnGeneMat(PanGEM_RNidx,:));
                temp_grRules5(PanGEM_RNidx) = strcat(temp_grRules5(PanGEM_RNidx),';',KO_stillNoHits{a});
            end
        end
    end
end

% Wow, basically no hits... just one corresponding with R00275
% (superoxide:superoxide oxidoreductase). I guess that means all the new
% KO's have missing reactions in PanGEM... yikes

% convenient table for doing manual inspection of these noHits KO's and
% their associated annotations, EC numbers and RNs.
ManualInspection = table(col1',col2',col3',col4');

% put together a list of KO's to keep after manual inspection. 
keepers = readtable(FileNames.keepers,'Delimiter',',','ReadVariableNames',false);

% ok now let's try to match this subset back to reactions, using col4,
% PanGEM.rnxs, and manual matching where necessary

% keepers subset first
for a = 1:numel(keepers.Var1)
    tempKO = keepers.Var1{a};
    colIdx = find(strcmp(tempKO,KO_stillNoHits));
    keepers_RN{a} = col4{colIdx};
end

% now the subset that had hits to KEGG
hitsIdx = find(nHits);
KO_hits = Pro_Assembly.uniqueKO(hitsIdx);
for a = 1:numel(hitsIdx)
    RN_hits{a} = RN_all{hitsIdx(a)};
end

% and finally the subset that had no hits to KEGG but had hits to PanGEM

for a = 1:numel(KO_noHits_keep)
    tempKO = KO_noHits_keep{a};
    PanGEM_idx = find(strcmp(tempKO,PanGEM.genes));
    RN_noHits_keep{a} = PanGEM.rxns(find(PanGEM.rxnGeneMat(:,PanGEM_idx)));
end

% so now let's merge the all the keepers and their reactions
fullKO_set = [keepers.Var1; KO_hits; KO_noHits_keep]    
fullRN_set = [keepers_RN'; RN_hits'; RN_noHits_keep'];

% let's see how many of these are already in PanGEM, which are to be
% removed, and which are to be added
noChange_set = intersect(fullKO_set,PanGEM.genes);
remove_set = setdiff(PanGEM.genes, fullKO_set);
add_set = setdiff(fullKO_set,PanGEM.genes);

% Let's manually inspect remove_set, since it's short we should be able
% to see if removing any of these will be an issue... ok did that... here's
% a file with the ones that should be able to be safely removed:

toRemove = readtable(FileNames.KO_remove,'Delimiter',',','ReadVariableNames',true);
canRemove = toRemove.KO(find(toRemove.remove==1));
cannotRemove = toRemove.KO(find(toRemove.remove==0));

% This 'cannotRemove' subset is small (currently 19 KO's) so we'll deal
% with them one-by-one as we move along.

% From the canRemove subset, let's loop through and see if that means
% removing reactions and metabolites as well. Each time we should check
% that the model can still grow. 
tempMod = PanGEM;
for a = 1:numel(canRemove)
    tempKO = canRemove{a};
    KO_idx = find(strcmp(tempKO,tempMod.genes));
    tempRN_idx = find(tempMod.rxnGeneMat(:,KO_idx));
    tempRN = tempMod.rxns(tempRN_idx);
    possibleRemoveRN{a} = tempRN;
    if ~isempty(tempRN_idx)
        nRN = numel(tempRN_idx);
        for b = 1:nRN
            nGenes = numel(find(PanGEM.rxnGeneMat(find(strcmp(tempRN{b},tempMod.rxns)),:)));
            if nGenes ==1
                tempMod2 = tempMod;
                tempMod2.lb(find(strcmp(tempRN{b},tempMod.rxns))) = 0;
                tempMod2.ub(find(strcmp(tempRN{b},tempMod.rxns))) = 0;
                tempSol = solveLP(tempMod2);
                if tempSol.f < -0.075
                    okToRemove{a}(b) = 1;
                    tempMod = removeRxns(tempMod2,tempRN{b},true,false,false);
                else
                    okToRemove{a}(b) = 0;
                end
            else
                okToRemove{a}(b)= 1;
            end
        end
    else
        okToRemove{a} = 1;
    end
end

% now let's get rid of reactions, genes, and metabolites associated with
% okToRemove genes
tempMod = PanGEM;
for a =1:numel(okToRemove)
    nRN = numel(possibleRemoveRN{a});
    if sum(okToRemove{a}) == nRN
        okToRemove2(a) = 1;
    end
end

% ok let's remove them
[rmSet,junk2,rmSet_idx] = intersect(canRemove(find(okToRemove2)),PanGEM.genes);
PanGEM.genes(rmSet_idx) = [];

% so the remaining list of genes which couldn't be removed is:
couldntRemove = [cannotRemove; canRemove(find(~okToRemove2'))]



% now let's generate a grRules vector for all the genes we want to add
% back. To do that, let's generate a list of all the reactions associated
% with the noChange_set and the add_set, then go back through and link
% those up to their respective KO's. 

newKO_set = setdiff([noChange_set;add_set;couldntRemove],rmSet);
n_newKO_set = numel(newKO_set);
% get all unique reactions associated with newKO_set
for a = 1:n_newKO_set
    dbIdx = find(strcmp(newKO_set{a},KO_RN.KO));
    if ~isempty(dbIdx)
        n_dbIdx = numel(dbIdx);
        for b = 1:n_dbIdx
            newRN_set{a}(b) = KO_RN.RN(dbIdx(b));
        end
    end
end
unique_newRN_set = unique([newRN_set{:}]);

% now link those up to their KO's
for a = 1:numel(unique_newRN_set)
    dbKO = KO_RN.KO(find(strcmp(unique_newRN_set{a},KO_RN.RN)));
    RN_KO = intersect(dbKO,newKO_set);
    new_grRules{a} = strjoin(RN_KO,';');
end

% some of these are already in PanGEM, so let's just replace those
for a = 1:numel(PanGEM.rxns)
    tempRN = PanGEM.rxns{a};
    tempRN_idx = find(strcmp(tempRN,unique_newRN_set));
    if ~isempty(tempRN_idx)
        new_grRules_PanGEM{a} = new_grRules{tempRN_idx};
    else
        new_grRules_PanGEM{a} = PanGEM.grRules{a};
    end
end
% revert format excel
new_grRules_PanGEM = strrep(new_grRules_PanGEM,'(','');
new_grRules_PanGEM = strrep(new_grRules_PanGEM,')','');
new_grRules_PanGEM = strrep(new_grRules_PanGEM,' or ',';');

% the rest of unique_newRN_set will be candidate reactions to add
[addSet_RN, addSet_RN_idx] = setdiff(unique_newRN_set,PanGEM.rxns);
addSet_grRules = new_grRules(addSet_RN_idx);

% let's see how many KO's are left out if we don't add these reactions in
leftOutKOs = unique(strsplit(strjoin(addSet_grRules,';'),';'));
% lots (401) are left out... ok so I guess we should leaf through these...
% yikes.

% it will help to know how many strains have these candidate ko's
for a = 1:numel(addSet_grRules)
    tempKOs = strsplit(addSet_grRules{a},';');
    n_tempKOs = numel(tempKOs);
    for b = 1:n_tempKOs
        nStrains_KO(b) = numel(find(Pro_Assembly.PresenceAbsenceMatrix(find(strcmp(tempKOs{b},Pro_Assembly.uniqueKO)),:)));
    end
    nStrains_KO2(a) = max(nStrains_KO);
    clear nStrains_KO
end



% let's at least get the equations for these
keggurl = 'http://rest.kegg.jp/get/'
clear response
for i = 1:numel(addSet_RN)
    extension = [keggurl,addSet_RN{i}];
    response{i} = urlread(extension);
end

for i = 1:numel(response)
    parsedResponse = textscan(response{i},'%s','delimiter','\t');
    nsearch = numel(parsedResponse{1});
    % get index for name, definition, equation, enzyme (EC), pathway
    nameInd = strfind(parsedResponse{1},'NAME');
    nameInd2 = find(~cellfun(@isempty,nameInd));
    definitionInd = strfind(parsedResponse{1},'DEFINITION');
    definitionInd2 = find(~cellfun(@isempty,definitionInd));
    equationInd = strfind(parsedResponse{1},'EQUATION');
    equationInd2 = find(~cellfun(@isempty,equationInd));
    enzymeInd = strfind(parsedResponse{1},'ENZYME');
    enzymeInd2 = find(~cellfun(@isempty,enzymeInd));
    pathwayInd = strfind(parsedResponse{1},'PATHWAY');
    pathwayInd2 = find(~cellfun(@isempty,pathwayInd));
    % assign to rxnStruct
    rxnStruct(i).NAME = strtrim(strrep(parsedResponse{1}(nameInd2),'NAME',''));
    rxnStruct(i).DEFINITION = strtrim(strrep(parsedResponse{1}(definitionInd2),'DEFINITION',''));
    rxnStruct(i).EQUATION = strtrim(strrep(parsedResponse{1}(equationInd2),'EQUATION',''));
    rxnStruct(i).ENZYME = strtrim(strrep(parsedResponse{1}(enzymeInd2),'ENZYME',''));
    rxnStruct(i).PATHWAY = strtrim(strrep(parsedResponse{1}(pathwayInd2),'PATHWAY',''));
end


for a = 1:numel(rxnStruct)
    equation = cell2mat(rxnStruct(a).EQUATION);
    metInd = regexp(equation,'C');
    if ~isempty(metInd)
        for j = 1:numel(metInd)
            metID{j} = equation(metInd(j):metInd(j)+5);
        end
        metID2{a} = metID;
    else metID2{a} = {};
        
    end
    clear metID
end

response = urlread('http://rest.kegg.jp/list/cpd');
response2 = textscan(response,'%s%s', 'Delimiter','\t');
cpdDat.Var1 = strrep(response2{1},'cpd:','');
cpdDat.Var2 = response2{2};

% The first metabolite name listed in KEGG is used, these are not IUPAC
% names, but rather many abbreviations and common names. Symbol
% substitutions to make importing with importExcelModel.m in RAVEN easier
% are as follows:

% ' ' -> '_'
% ',' -> ''
% ''' -> ''
% '-' -> '_'
% '(' -> ''
% ')' -> ''
% '+' -> ''
% '[' -> ''
% ']' -> ''
% '/' -> '_'
% '>' -> ''
% '"' -> ''
% 'acyl_carrier_protein' -> 'acp'

% after these substitutions, we remove double underscores and replace with
% single
% '__' -> '_'

for a = 1:numel(cpdDat.Var2)
    cpdStr = strsplit(cpdDat.Var2{a},';');
    cpdStr2 = cpdStr{1};
    cpdStr3 = strrep(cpdStr2,' ','_');
    cpdStr4 = strrep(cpdStr3,',','');
    cpdStr5 = strrep(cpdStr4,'-','_');
    cpdStr6 = strrep(cpdStr5,'(','');
    cpdStr7 = strrep(cpdStr6,')','');
    cpdStr8 = strrep(cpdStr7,'+','');
    cpdStr9 = strrep(cpdStr8,'[','');
    cpdStr10 = strrep(cpdStr9,']','');
    cpdStr11 = strrep(cpdStr10,'/','');
    cpdStr12 = strrep(cpdStr11,'>','');
    cpdStr13 = strrep(cpdStr12,'"','');
    cpdStr14 = strrep(cpdStr13,'acyl_carrier_protein','acp');
    cpdStr15 = strrep(cpdStr14,'__','_');
    cpdDat.Var2{a} = cpdStr15;
end

% replace cpd ID's with new names
for a = 1:numel(rxnStruct)
    EQmets = metID2{a};
    equation = cell2mat(rxnStruct(a).EQUATION);
    metInd = regexp(equation,'C');
    for j = 1:numel(EQmets)
        cpdNameInd(j) = find(strcmp(EQmets{j},cpdDat.Var1));
        equation = strrep(equation,EQmets{j},cpdDat.Var2{cpdNameInd(j)});
    end
    newEquation{a} = equation;
end

% update addSet structure
addSet = struct2table(rxnStruct);
addSet.RN = addSet_RN';
addSet.grRules = addSet_grRules';
addSet.EQUATIONverbose = newEquation';
addSet.nStrains_grRules = nStrains_KO2';

% we've whittled down the list of new reactions to something manageable.
% Now let's loop through, adding a new reaction each time and seeing if it
% causes loops or other changes. If not, we can keep it, if so, we have to
% look for directionality corrections.

% I've gone through an whittled this subset down even further by manual
% inspection. So let's find the indices of these:
fileName = 'Users/jrcasey/Documents/MATLAB/CBIOMES/Pangenomes/Prochlorococcus/toAdd_20200823.csv';
toAdd = readtable(fileName,'Delimiter',',','ReadVariableNames',true);
for a = 1:numel(toAdd.RN)
    addSet_idx(a) = find(strcmp(toAdd.RN{a},addSet.RN));
end


% first let's go through one-by-one

tempMod = PanGEM;
for a = 1:numel(toAdd.RN)
    rxnsToAdd.rxns = {toAdd.RN{a}};
    rxnsToAdd.equations = {addSet.EQUATIONverbose{addSet_idx(a)}};
    rxnsToAdd.lb = toAdd.LB(a);
    rxnsToAdd.ub = toAdd.UB(a);
    tempMod2 = addRxns(tempMod,rxnsToAdd,1,'c',true);
    tempSol = solveLP(tempMod2,1);
    if tempSol.stat
        mu(a) = -tempSol.f;
    end
end






%% Generate Gene-reaction matrix for PanGEM
grMat = zeros(numel(PanGEM.rxns),numel(Pro_Assembly.uniqueKO));
for a = 1:numel(Pro_Assembly.uniqueKO)
    tempRN = RN_all{a};
    if ~isempty(tempRN)
        nRN = numel(tempRN);
        for b = 1:nRN
            RN_idx = find(strcmp(tempRN{b},PanGEM.rxns));
            grMat(RN_idx,a) = 1;
        end
    end
end

    



%% Get reactions to remove from PanGEM
removeRN = setdiff(PanGEM.rxns,uniqueRN);
% lets see which of those cannot carry flux anyway
removeRN_flux = haveFlux(PanGEM,1e-12,removeRN)



% from this list, check whether another KO which is found in the
% Pro_Assembly database can be associated with that reaction
RN = 'R00028';
tempResponse = urlread(strcat('http://rest.kegg.jp/get/',RN));


% check whether an associated KO is found in the Pro_Assembly and in the
% KO_RN database
temp = 'K01187';
find(strcmp(temp,Pro_Assembly.uniqueKO)) 
find(strcmp(temp,KO_RN.KO))
tempResponse = urlread(strcat('http://rest.kegg.jp/get/',temp));
% just print the response up until the GENES line
GENES_idx = regexp(tempResponse,'GENES');
tempResponse(1:GENES_idx)


%% Map rxns to strains
% for a = 1:nStr
%     strKO = Pro_Assembly.Strains.(strains{a}).KO;
%     % find indices of strain KO's in uniqueKO
%     [junk,strKO_idx, uniqueKO_idx] = intersect(strKO,Pro_Assembly.uniqueKO);
%     
%     Pro_Assembly.Strains.(strains{a}).RN = 0;

%% Map new KO's to PanGEM rxns... see what's left over


%% Matrix of rxns vs strains
% preallocate
PAMat_RN = zeros(n_uniqueRN,nStr);
for a = 1:nStr
    strRN = Pro_Assembly.Strains.(strains{a}).KO;
    for b = 1:n_uniqueRN
        PAMat_RN(b,a) = ismember(uniqueRN{b},strKO);
    end
end

