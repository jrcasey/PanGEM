%% Step 8: proteomeGEM (optional)
% Imports a GEM generated using PanGEM Toolbox, imports the corresponding
% proteome sequence, and generates synthesis reactions for each enzyme.
% Saves the appended GEM as a mat file for implementation with MSE.


%% Directories
FileNames.orgDB = 'data/genomes/orgDatabase_20200805.csv';
FileNames.strainModelsDir = 'data/models/strains/';
FileNames.Assembly_L1 = 'data/assemblies/Assembly_L1_20200830.mat';

%% Import data
% Load L1 assembly
load(FileNames.Assembly_L1) % variable name is Pro_Assembly_L1

% import orgDB
orgDB = readtable(FileNames.orgDB,'Delimiter',',','ReadVariableNames',true);

% get list of strains to import
strList = dir(FileNames.strainModelsDir);
strList = {strList.name};
strList([1 2]) = []; % remove . and .. directories
if any(contains(strList,'.DS_Store'))
    strList(find(contains(strList,'.DS_Store'))) = [];
end
strList = strrep(strList,'.mat','');
nStr = numel(strList);


%% Compile proteomeGEM for each strain
% loop through strains, loading each strain GEM and adding synthesis
% reactions for each enzym and the corresponding enzymes as metabolites.
for a = 1:nStr
    strName = strList{a};
    % load strMod
    load(strcat(FileNames.strainModelsDir,strName,'.mat'));
    model = addEnzymeSynthesis(model,Pro_Assembly_L1);
    save(['data/models/pStrains/pGEM_',strName,'.mat'],'model');
end

%% Solve optimization

% k_cat for each enzyme
k_cat_conv = (1/6.022e23) .* 1000 .* 3600; % conversion factor from molecules enzyme-1 s-1 to mmol enzyme-1 h-1
k_cat_c = 10; % molecules enzyme-1  s-1, a constant as a guess until we have these data
k_cat = k_cat_conv .* repmat(k_cat_c,numel(model.genes),1); % mmol enzyme-1 h-1


% Compute initial point
tempSol = solveLP(model,1); 
% convert to g of enzyme per gDW
for a = 1:numel(model.rxns)
    rxnGenes_idx = find(model.rxnGeneMat(a,:));
    % just pick one for now
    if ~isempty(rxnGenes_idx)
        for b = 1:numel(rxnGenes_idx)
            E_gDW(b) = abs(tempSol.x(a)) ./ k_cat(rxnGenes_idx(1)); % enzymes gDW-1
        end
        E_g_gDW(a) = nansum(E_gDW' .* model.geneProductMW(rxnGenes_idx)) .* (1/6.022e23);
    else
        E_g_gDW(a) = NaN;
    end
    clear E_gDW
end

total_proteome = nansum(E_g_gDW);

% compute total for each enzyme
for a = 1:numel(model.genes)
    geneRxns_idx = find(model.rxnGeneMat(:,a));
    total_enzyme(a) = nansum(E_g_gDW(geneRxns_idx));
end

% % top50
% [total_enzyme_sorted, sortOrder] = sort(total_enzyme,'descend')
% table(model.genes(sortOrder(1:50)),total_enzyme_sorted(1:50)')

% assign to initial vector
x0 = total_enzyme;

lb = zeros(numel(model.genes),1);
ub = repmat(0.1,numel(model.genes),1);

A = ones(numel(model.genes),1);
b = 0.5; % maximum proteome g gDW-1

% set up optimization problem
options = optimoptions('fmincon','ConstraintTolerance',1e-3,'MaxIterations',1000);
prob = struct;
prob.objective = @(x)solvePFBA(x,model,k_cat);
prob.Aineq = A;
prob.bineq = b;
prob.Aeq = [];
prob.beq = [];
prob.lb = lb;
prob.ub = ub;
prob.nonlcon = [];
prob.solver = 'fmincon';
prob.options = options;


