%% Step 4 - Annotate manually curated excel model
% This step loads an excel model (see ref) and converts it to the more
% widely adopted COBRA format. The COBRA model is then annotated using a
% variety of common SysBio databases for easier sharing.

% (About 30 minutes runtime)

%% Generate COBRA model from Excel model
excelModPath = 'data/models/PanGEM_pro_20200817.xlsx';
PanGEM = excelToCobra(excelModPath);

% check that the model is feasible
changeCobraSolver('ibm_cplex');
sol = optimizeCbModel(PanGEM,[],'one');

%% Annotate model
LinkDB_Dir = 'data/dbLinks/';
PanGEM = annotateModel(PanGEM,LinkDB_Dir);

%% Save annotated model to mat and/or SBML formats
PanGEM.description = 'Step 4 PanGEM'
% mat file (quickest)
save('data/models/Pro_PanGEM.mat','PanGEM');

% sbml (for sharing)
% SBML_Model = writeCbModel(PanGEM,'format','sbml','fileName','data/models/Pro_PanGEM.sbml');
