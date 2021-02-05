function cobraMod = excelToCobra(excelModPath,cobraModPath)
%% Convert an excel model to COBRA using RAVEN
% 
% Inputs
% excelModPath            String. Path to excel model
% cobraModPath            String. Path to cobra SBML file output (optional)
% 
% Outputs
% cobraMod                Structure. Cobra model.
% 
% John R. Casey 
% 08/26/2020

%% Check if cobraModPath is specified
if nargin == 2
    saveMod = 1;
else
    saveMod = 0;
end

%% import excel model with RAVEN
fileName = sprintf('%s',excelModPath);
ravenMod = importExcelModel(fileName,false,true,false);

% Remove exchange reactions for COBRA
ravenMod2 = simplifyModel(ravenMod);

% save raven mat file to excelModPath folder
[filepath,name,ext] = fileparts(excelModPath)
ravenModPath = strcat(filepath,'/',name,'_raven.mat');
save(ravenModPath,'ravenMod2');

%% import with COBRA
cobraMod = readCbModel(ravenModPath);


%% Save SBML file if prompted
if saveMod
    SBML_Model = writeCbModel(cobraMod,'format','sbml','fileName',cobraModPath)
end

end

