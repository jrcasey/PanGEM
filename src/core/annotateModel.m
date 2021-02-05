function annotatedModel = annotateModel(cobraMod,LinkDB_Dir)

%% Annotate Reactions
% add KEGG and MetaNetX ID for reactions

% add KEGG reaction ID's to model
url = 'http://rest.genome.jp/list/rn';
response = urlread(url);
response2 = textscan(response,'%s %s','Delimiter','\t');
rnList = strrep(response2{1},'rn:','');
for a = 1:numel(cobraMod.rxns)
    match = ismember(cobraMod.rxns{a},rnList);
    if match
        cobraMod.rxnKEGGID{a} = cobraMod.rxns{a};
    end
end

% add MetaNetX ID's to model
fileName = strcat(LinkDB_Dir,'reac_xref.txt');
MetaNetX_rxn = readtable(fileName,'Delimiter','\t','ReadVariableNames',false,'Format','%s%s%s');
cobraMod.rxnMetaNetXID = cell(numel(cobraMod.rxns),1);
for a = 1:numel(cobraMod.rxns)
    tempKEGG = cobraMod.rxnKEGGID{a};
    if ~isempty(tempKEGG)
        newKEGG = strcat('kegg:',tempKEGG);
        datIdx = find(strcmp(newKEGG,MetaNetX_rxn.Var1));
        if ~isempty(datIdx)
            cobraMod.rxnMetaNetXID{a} = MetaNetX_rxn.Var2{datIdx(1)};
        end
    end
end
clear MetaNetX_rxn

%% Annotate Metabolites
% Add PubChem, Chebi, InChI, MetaNetX IDs to metabolites

% start with a list of KEGG cpd ID's from the model. You'll need to geneate
% a separate csv file for this task. Note that you need to trim the
% boundary metabolites from this list. 
fileName = strcat(LinkDB_Dir,'PanGEM_pro_met_KEGG.csv');
PanGEM_metKEGG = readtable(fileName,'Delimiter',',','ReadVariableNames',true);
% get indices of this database to the model
for a = 1:numel(PanGEM_metKEGG.ID)
    PanGEM_metIdx(a) = find(strcmp(PanGEM_metKEGG.ID{a},cobraMod.mets));
end
% assign to model
cobraMod.metKEGGID = cell(numel(cobraMod.mets),1);
cobraMod.metKEGGID(PanGEM_metIdx) = PanGEM_metKEGG.KEGG_cpd;

% Add PubChem ID's
url = 'http://rest.genome.jp/link/cpd/pubchem';
response = urlread(url);
response2 = textscan(response,'%s %s %s','Delimiter','\t');
cpdList = strrep(response2{2},'cpd:','');
pubChemList = strrep(response2{1},'pubchem:','');

cobraMod.metPubChemID = cell(numel(cobraMod.mets),1);
for a = 1:numel(PanGEM_metKEGG.KEGG_cpd)
    tempcpd = PanGEM_metKEGG.KEGG_cpd{a};
    if ~isempty(tempcpd)
        dbIdx = find(strcmp(PanGEM_metKEGG.KEGG_cpd{a},cpdList));
        cobraMod.metPubChemID(PanGEM_metIdx(a)) = pubChemList(dbIdx);
    end
end

% add CheBI ID's to mets
fileName = strcat(LinkDB_Dir,'compound_chebi.csv');
ChEBI_KEGG = readtable(fileName,'Delimiter',',','ReadVariableNames',false);

cobraMod.metChEBIID = cell(numel(cobraMod.mets),1);
for a = 1:numel(cobraMod.mets)
    dbIdx = find(strcmp(cobraMod.metKEGGID{a},ChEBI_KEGG.Var1));
    if ~isempty(dbIdx)
        cobraMod.metChEBIID{a} = ChEBI_KEGG.Var2(dbIdx(1));
    end
end

% add InChI strings for each metabolite
cobraMod.metInChI = cell(numel(cobraMod.mets),1);
urlPrefix = 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:';
for a = 1:numel(cobraMod.mets)
    ChEBI_ID = mat2str(cobraMod.metChEBIID{a});
    if ~isempty(cobraMod.metChEBIID{a})
        temp = urlread(strcat(urlPrefix,ChEBI_ID));
        startIdx = regexp(temp,'<td>InChI=');
        if ~isempty(startIdx)
            allEndIdx1 = regexp(temp,'</td>');

            % only take the indices larger than startIdx
            allEndIdx1(find(allEndIdx1<startIdx)) = [];
           
            endIdx = allEndIdx1(1);
            tempInChI = temp(startIdx+4:endIdx-1);
            tempInChI2 = erase(tempInChI,'<img alt="" src="images/invisible_space.gif" width="0" height="0" border="0"/>');
            cobraMod.metInChI{a} = tempInChI2;
        end
    end    
end

% add MetaNetX IDs
fileName = strcat(LinkDB_Dir,'chem_xref.txt');
MetaNetX_cpd = readtable(fileName,'Delimiter','\t','ReadVariableNames',false,'Format','%s%s%s');
cobraMod.metMetaNetXID = cell(numel(cobraMod.mets),1);
for a = 1:numel(cobraMod.mets)
    tempChEBI = cobraMod.metChEBIID{a};
    if ~isempty(tempChEBI)
        newChEBI = strcat('chebi:',mat2str(tempChEBI));
        datIdx = find(strcmp(newChEBI,MetaNetX_cpd.Var1));
        if ~isempty(datIdx)
            cobraMod.metMetaNetXID{a} = MetaNetX_cpd.Var2{datIdx(1)};
        end
    end
end
clear MetaNetX_cpd

%% Add SBO Terms to all model features

% Add SBO terms to metabolites
cobraMod.metSBOTerms = cell(numel(cobraMod.mets),1);
cobraMod.metSBOTerms(:) = {'SBO:0000247'};

% Add SBO terms to reactions
cobraMod.rxnSBOTerms = cell(numel(cobraMod.rxns),1);
cobraMod.rxnSBOTerms(:) = {'SBO:0000176'};
% need to assign a different value to exchange reactions
cobraMod.subSystems = cellfun(@cell2str,cobraMod.subSystems,'UniformOutput',false); % something weird with the subSystems array
ExIdx = find(strcmp('Exchange',cobraMod.subSystems));
cobraMod.rxnSBOTerms(ExIdx) = {'SBO:0000627'};
% and another to transport reactions
ExIdx = find(strcmp('Transport',cobraMod.subSystems));
cobraMod.rxnSBOTerms(ExIdx) = {'SBO:0000185'};
% finally one for the BOF
bofIdx = find(cobraMod.c);
cobraMod.rxnSBOTerms(bofIdx) = {'SBO:0000629'};

% Add SBO terms to genes
cobraMod.geneSBOTerms = cell(numel(cobraMod.genes),1);
cobraMod.geneSBOTerms(:) = {'SBO:0000243'};

%% Store output
annotatedModel = cobraMod;
annotatedModel.modelID = 'COBRA Model - PanGEM comprised of 649 Prochlorococcus genomes'

end
