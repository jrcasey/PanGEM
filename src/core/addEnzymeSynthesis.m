function model = addEnzymeSynthesis(model,Assembly_L1)
% This function searches the level 1 assembly for genes corresponding to
% each KEGG Ortholog present in a strain GEM. It then adds the
% corresponding gene product 
strName = model.id;
nGenes = numel(model.genes);
model.geneProductMW = NaN(nGenes,1);
for b = 1:nGenes
    % AA sequence for each KO
    KO = model.genes{b};
    % find corresponding gene(s)
    KOidx = find(strcmp(KO,Assembly_L1.Strains.(strName).KO));
    if ~isempty(KOidx)
        AAseq = Assembly_L1.Strains.(strName).Sequence{KOidx(1)};
        AAprofile = getAADistribution(AAseq); % get amino acid composition of protein
        AAstats = getProteinElemental(AAprofile); % get MW and elemental composition of protein
        model.geneProductMW(b) = AAstats.MW;
        % construct equation
        tmpStr = '';
        for c = 1:numel(AAprofile.names)
            tmpStr = [tmpStr,' + ',num2str(1e-3.*AAprofile.count(c)),' ',AAprofile.names{c}]; % try scaling down the coefficients to avoid fluxes less than the solver tolerance
        end
        tmpStr([1:3]) = [];
        tmpStr = [tmpStr,' => 0.001 E_',KO];
        
        % add enzyme as 'metabolite'
        metsToAdd = struct;
        metsToAdd.mets = {strcat('E_',KO)};
        metsToAdd.metNames = {['Enzyme corresponding to',' ', KO,' ','based on sequence from',' ',Assembly_L1.Strains.(strName).GeneID{KOidx(1)}]};
        metsToAdd.compartments = {'c'};
        metsToAdd.formula = {AAstats.formula};
        model = addMets(model,metsToAdd);
        model.metCharges(end+1) = NaN;
        model.metKEGGID{end+1} = {};
        model.metChEBIID{end+1} = {};
        model.metInChI{end+1} = {};
        model.metMetaNetXID{end+1} = {};
        model.metPubChemID{end+1} = {};
        
        % add enzyme synthesis reaction
        rxnsToAdd = struct;
        rxnsToAdd.rxns = {['E_',KO,'_synthesis']};
        rxnsToAdd.rxnNames = {['Synthesis of enzyme corresponding to',' ', KO,' ','based on sequence from',' ',Assembly_L1.Strains.(strName).GeneID{KOidx(1)}]};
        rxnsToAdd.equations = {tmpStr};
        rxnsToAdd.lb = 0;
        rxnsToAdd.ub = 1000;
        rxnsToAdd.subSystems = {'Enzyme Synthesis'};
        model = addRxns(model,rxnsToAdd,1,'c');
        model.orthologPresent(b) = 1;
        
        % add gene product elemental stats
        model.geneProduct_NitrogenAtoms(b) = AAstats.N;
        model.geneProduct_CN(b) = AAstats.CN;
    else
        model.orthologPresent(b) = 0;
        
        model.geneProduct_NitrogenAtoms(b) = NaN;
        model.geneProduct_CN(b) = NaN;
    end
end


end
