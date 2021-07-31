function ProteinElemental = getProteinElemental(AA_Distribution)

% Retrieves the molecular weight for a protein.

% Inputs
% AA_Distribution - Structure. Output from getAADistribution.

% Outputs
% LocusMW         - Double. Molecular weight

% Get elemental composition and molecular weight for each amino acid (be
% sure this is in the same order as AA_Distribution.names
AAFormulas = {'C3H7NO2';'C6H14N4O2';'C4H8N2O3';'C4H7NO4';'C3H7NO2S';'C5H9NO4';...
    'C5H10N2O3';'C2H5NO2';'C6H9N3O2';'C6H13NO2';'C6H13NO2';'C6H14N2O2';'C5H11NO2S';...
    'C9H11NO2';'C5H9NO2';'C3H7NO3';'C4H9NO3';'C11H12N2O2';'C9H11NO3';'C5H11NO2'};

for i = 1:numel(AAFormulas)
    [elements, useMat{i}, exitFlag{i}, MW{i}] = parseFormulas(AAFormulas(i), true,false,true);
end

% protein molecular formula
useMat2 = cell2mat(useMat');
prtFormula_dat = useMat2.*repmat(AA_Distribution.count',1,numel(elements.abbrevs));
prtFormula_Sum = sum(prtFormula_dat,1);
% get element indices
Cind = find(strcmp('C',elements.abbrevs));
Hind = find(strcmp('H',elements.abbrevs));
Nind = find(strcmp('N',elements.abbrevs));
Oind = find(strcmp('O',elements.abbrevs));
Pind = find(strcmp('P',elements.abbrevs));
Sind = find(strcmp('S',elements.abbrevs));
Seind = find(strcmp('Se',elements.abbrevs));

ProteinElemental.formula = strcat('C',mat2str(prtFormula_Sum(Cind)),'H',mat2str(prtFormula_Sum(Hind)),...
    'N',mat2str(prtFormula_Sum(Nind)),'O',mat2str(prtFormula_Sum(Oind)),...
    'P',mat2str(prtFormula_Sum(Pind)),'S',mat2str(prtFormula_Sum(Sind)),...
    'Se',mat2str(prtFormula_Sum(Seind)));

MW2 = [MW{:}];
ProteinElemental.MW = sum(MW2.*AA_Distribution.count);
ProteinElemental.CN = prtFormula_Sum(Cind)./prtFormula_Sum(Nind);
ProteinElemental.N = prtFormula_Sum(Nind);
end

