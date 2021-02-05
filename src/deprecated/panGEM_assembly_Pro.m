%% panGEM assembly
% This script takes the genes and KO's generated in Pangenome_assembly_redo
% and populates a GEM model structure. KO's are first linked to KEGG
% reactions (R#####) which are in turn linked to metabolites. Genes and strains from
% the database are then added to the structure as grRules and Strains.
% We'll later use Strains to generate strain GEMs but for now, let's create
% a panGEM!

% First let's import the database
load('Pro_Assembly2.mat');

%% Import reaction associations to KO's from KEGG REST API relational database
% use linkDB database for matching, that way we don't have to
% urlread every entry. 
temp = strrep(Dat.KO,'NaN','');
Dat.KO = temp;
clear temp;

% import relational database
databasepath = 'CBIOMES/Scripts/KEGG_KO_RN.csv';
KO_RN = readtable(databasepath,'delimiter',',','ReadVariableNames',false);

for i = 1:numel(Dat.KO);
    if ~isempty(Dat.KO{i});
        matchInd = find(strcmp(Dat.KO{i},KO_RN.Var1));
        matchIndn = numel(matchInd);
        listRxn{i} = KO_RN.Var2(matchInd);
        listRxnInd{i} = repmat({i},matchIndn,1);
    else listRxn{i} = [];
    end
end

Dat.RXNS = listRxn';
Dat.RXNSind = listRxnInd';
% save and clear to free up memory
savefast('CBIOMES/Pangenomes/Prochlorococcus/Pro_Assembly3.mat','Dat')
clear

%% Link Reactions to equations

% use keggRxns.mat database, but first check how many reactions in our
% database match to it.

load('Pro_Assembly3.mat')
load('keggRxns.mat')

% alright now let's start an excel model, but not from scratch thank
% christ, import the MED4 model and see how much we get. 

fileName = sprintf('%s','Chalmers/Datasets/MED4_models/draftMED4/editing/draftMED4_20161114Combined.xlsx')
model = importExcelModel(fileName,true,true,false);
unused = model.mets(~any(model.S,2));
model = removeMets(model,unused);


vertRxns = vertcat(Dat.RXNS{:});
vertRxnInd = vertcat(Dat.RXNSind{:});
uniqueRxns = unique(vertRxns);

overlap = intersect(uniqueRxns,model.rxns);
missing = setdiff(uniqueRxns,model.rxns);
% not horrible, but not great either. Will have to add these manually
% somehow.

%% TIME TO GO MANUAL! LUL! ughhhhhhhh

% first we need to set up as many columns as possible. The UID is RXNS
% number, so let's start from Dat and get a list of all KO's, genes, and
% strains associated with each R number. We can use unique RXNS to get that
% going. 
rxnStruct = struct;
for i = 1:numel(uniqueRxns);
    rxnInd = find(strcmp(uniqueRxns{i},vertRxns));
    rxnStruct.GeneID{i} = strjoin(Dat.GeneID(cell2mat(vertRxnInd(rxnInd)))',';')';
    rxnStruct.Strain{i} = strjoin(Dat.Strain(cell2mat(vertRxnInd(rxnInd)))',';')';
    rxnStruct.KO{i} = strjoin(Dat.KO(cell2mat(vertRxnInd(rxnInd)))',';')';
end


for i = 1:numel(model.rxns);
    keep(i) = any(strcmp(model.rxns{i},uniqueRxns));
end
keep = keep';


%% ugh nah, let's retrieve the entries from the KEGG API
keggurl = 'http://rest.kegg.jp/get/'
for i = 1:numel(uniqueRxns);
    extension = [keggurl,uniqueRxns{i}];
    response{i} = urlread(extension);
end

for i = 1:numel(response);
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
    rxnStruct.NAME{i} = strtrim(strrep(parsedResponse{1}(nameInd2),'NAME',''));
    rxnStruct.DEFINITION{i} = strtrim(strrep(parsedResponse{1}(definitionInd2),'DEFINITION',''));
    rxnStruct.EQUATION{i} = strtrim(strrep(parsedResponse{1}(equationInd2),'EQUATION',''));
    rxnStruct.ENZYME{i} = strtrim(strrep(parsedResponse{1}(enzymeInd2),'ENZYME',''));
    rxnStruct.PATHWAY{i} = strtrim(strrep(parsedResponse{1}(pathwayInd2),'PATHWAY',''));
end

rxnStruct.GeneID = rxnStruct.GeneID';
rxnStruct.Strain = rxnStruct.Strain';
rxnStruct.KO = rxnStruct.KO';

rxnStruct.NAME = rxnStruct.NAME';
rxnStruct.DEFINITION = rxnStruct.DEFINITION';
rxnStruct.EQUATION = rxnStruct.EQUATION';
rxnStruct.ENZYME = rxnStruct.ENZYME';
rxnStruct.PATHWAY = rxnStruct.PATHWAY';
rxnStruct.RXNS = uniqueRxns;

% add lb, ub, and compartment

for i = 1:numel(uniqueRxns);
    modInd = find(strcmp(uniqueRxns{i},model.rxns));
    if isempty(modInd);
        modInd2(i) = 2;
    else
        modInd2(i) = modInd;
end
end
rxnStruct.lb = model.lb(modInd2);
rxnStruct.ub = model.ub(modInd2);
rxnStruct.rxnComps = model.rxnComps(modInd2);

%% get metabolite names from compound ID's in reaction equations in rxnStruct
% also need to edit the names to be able to import with RAVEN

for i = 1:numel(rxnStruct.EQUATION);
    equation = cell2mat(rxnStruct.EQUATION{i});
    metInd = regexp(equation,'C');
    if ~isempty(metInd);
        for j = 1:numel(metInd);
            metID{j} = equation(metInd(j):metInd(j)+5);
        end
        metID2{i} = metID;
    else metID2{i} = [];
        
    end
    clear metID
end

% import relational table (retrieved from 'http://rest.kegg.jp/list/cpd')

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
% 'acyl_carrier_protein' -> 'acp'

% after these substitutions, we remove double underscores and replace with
% single
% '__' -> '_'

cpdPath = 'CBIOMES/Scripts/KEGG_cpd_list.csv'
cpdDat = readtable(cpdPath,'delimiter',',','ReadVariableNames',false);


% replace cpd ID's with new names
for i = 1:numel(rxnStruct.EQUATION);
    EQmets = metID2{i};
    equation = cell2mat(rxnStruct.EQUATION{i});
    metInd = regexp(equation,'C');
    for j = 1:numel(EQmets);
        cpdNameInd(j) = find(strcmp(EQmets{j},cpdDat.Var1));
        equation = strrep(equation,EQmets{j},cpdDat.Var2{cpdNameInd(j)});
    end
    newEquation{i} = equation;
end

% update rxnStruct
rxnStruct.EQUATIONverbose = newEquation';

% fill in metabolite sheet

% stack and unique the cpd equations
uniqueMets = unique([metID2{:}]);
for i = 1:numel(uniqueMets);
    tempind = find(strcmp(uniqueMets{i},cpdDat.Var1));
    uniqueMetsVerbose{i} = cpdDat.Var2{tempind};
end
metStruct = struct;
metStruct.metID = uniqueMets';
metStruct.metName = uniqueMetsVerbose';

% link cpd ID's to MED4 model fields

% first let's see about the coverage
metOverlap = intersect(model.mets,metStruct.metName)
metMissing = setdiff(model.mets,metStruct.metName)

for i = 1:numel(metStruct.metName);
    modMetInd = find(strcmp(metStruct.metName{i},model.mets));
    if ~isempty(modMetInd);
        metStruct.metFormulas{i} = model.metFormulas{modMetInd};
    else metStruct.metFormulas{i} = [];
    end
end

metStruct.metFormulas = metStruct.metFormulas';


% retrieve additional compound information from the KEGG API

keggurl = 'http://rest.kegg.jp/get/'

for i = 1:numel(metStruct.metID);
    extension = [keggurl,metStruct.metID{i}];
    metResponse{i} = urlread(extension);
end

for i = 1:numel(metResponse);
    parsedResponse = textscan(metResponse{i},'%s','delimiter','\t');
    
    % get index for FORMULA, MOL_WEIGHT, PubChem, and ChEBI
    formulaInd = strfind(parsedResponse{1},'FORMULA');
    formulaInd2 = find(~cellfun(@isempty,formulaInd));
    molwtInd = strfind(parsedResponse{1},'MOL_WEIGHT');
    molwtInd2 = find(~cellfun(@isempty,molwtInd));
    pubchemInd = strfind(parsedResponse{1},'PubChem');
    pubchemInd2 = find(~cellfun(@isempty,pubchemInd));
    chebiInd = strfind(parsedResponse{1},'ChEBI');
    chebiInd2 = find(~cellfun(@isempty,chebiInd));
    
    % assign to rxnStruct
    metStruct.FORMULA{i} = strtrim(strrep(parsedResponse{1}(formulaInd2),'FORMULA',''));
    metStruct.MOL_WEIGHT{i} = strtrim(strrep(parsedResponse{1}(molwtInd2),'MOL_WEIGHT',''));
    metStruct.PubChem{i} = strtrim(strrep(parsedResponse{1}(pubchemInd2),'PubChem',''));
    metStruct.ChEBI{i} = strtrim(strrep(parsedResponse{1}(chebiInd2),'ChEBI',''));

end
metStruct.FORMULA = metStruct.FORMULA';
metStruct.MOL_WEIGHT = metStruct.MOL_WEIGHT';
metStruct.PubChem = metStruct.PubChem';
metStruct.ChEBI = metStruct.ChEBI';

%% save rxnStruct
savefast('CBIOMES/Pangenomes/Prochlorococcus/rxnStruct.mat','rxnStruct')
savefast('CBIOMES/Pangenomes/Prochlorococcus/metStruct.mat','metStruct')

clearvars -except rxnStruct

%% start fresh
load('rxnStruct.mat');load('metStruct.mat');

%% Clean up rxnStruct.mat
for i = 1:numel(rxnStruct.GeneID);
    % Genes
    GeneID_split = strsplit(rxnStruct.GeneID{i}',';');
    GeneID_unique = unique(GeneID_split);
    newGeneID = strjoin(GeneID_unique,';');
    rxnStruct.GeneID{i} = newGeneID;
    clear GeneID_split GeneID_unique newGeneID;
    % Strains
    Strain_split = strsplit(rxnStruct.Strain{i}',';');
    Strain_unique = unique(Strain_split);
    newStrainID = strjoin(Strain_unique,';');
    rxnStruct.Strain{i} = newStrainID;
    clear Strain_split Strain_unique newStrainID;
    % KO's
    KO_split = strsplit(rxnStruct.KO{i}',';');
    KO_unique = unique(KO_split);
    newKOID = strjoin(KO_unique,';');
    rxnStruct.KO{i} = newKOID;
    clear Strain_split Strain_unique newStrainID;
end





%% Link RHEA reactions
RHEApath = 'CBIOMES/DatabaseLinks/RHEA/rhea2kegg_reaction.csv';
RHEAdat = readtable(RHEApath,'delimiter',',','ReadVariableNames',false,'Format','%f%f%s');
RHEAdat.Properties.VariableNames = {'RHEAID','MASTERID','KEGGReactionID'};

for i = 1:numel(rxnStruct.RXNS);
    RHEAind = find(strcmp(rxnStruct.RXNS{i},RHEAdat.KEGGReactionID));
    if ~isempty(RHEAind);
        RHEAID{i} = RHEAdat.RHEAID(RHEAind);
        MASTERID{i} = RHEAdat.MASTERID(RHEAind);
    else RHEAID{i} = [];
        MASTERID{i} = [];
    end
end
rxnStruct.RHEA = RHEAID'; % RHEA reaction ID (specified direction)
rxnStruct.RHEAMASTER = MASTERID'; % equivalent RHEA reaction ID (undefined direction)

%% Update rxnStruct
savefast('CBIOMES/Pangenomes/Prochlorococcus/rxnStruct.mat','rxnStruct')

%% Get Gibbs formation energy and charge from Jankowski et al., 2008
% Group contribution method

% Import data
dat_path = 'CBIOMES/Data/Jankowski_groupContributionMethod/dG_KEGG_CPD2.csv';
CPD_dG = readtable(dat_path,'Delimiter',',','ReadVariableNames',true)

% Match
for i = 1:numel(metStruct.metID)
    ind = find(strcmp(metStruct.metID{i},CPD_dG.ENTRY));
    temp_dG{i} = CPD_dG.DELTAG(ind);
    temp_uncertainty{i} = CPD_dG.UNCERTAINTY(ind);
    temp_charge{i} = CPD_dG.CHARGE(ind);
end

% update structure
metStruct.CHARGE = temp_charge';
metStruct.dG_formation = temp_dG';
metStruct.dG_formation_uncertainty = temp_uncertainty';

% Save
savefast('CBIOMES/Pangenomes/Prochlorococcus/metStruct.mat','metStruct')

%% Manual curation of rxnStruct - edited simultaneously in the excel file: 
% PanGEM_pro_201XXXXX.xlsx in CBIOMES/Pangenomes/Prochlorococcus/

% rxnStruct.var{ind} = 'value'
%ind = find






