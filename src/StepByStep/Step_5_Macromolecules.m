%% Step 5 - Macromolecules

% This step deals with the tricky business of macromolecular pools. SBML
% will only accept integer coefficients for metabolite formulas, so we have
% to do this separately.

% The main groups to attend to are acyl chains in the various lipids, and
% the macromolecular pools comprising biomass.

% (Less than 1 minute runtime)

%% Load annotated PanGEM
load('data/models/Pro_PanGEM.mat')

%% Acyl groups
% let's add acyl chain lengths to metabolites using the acyl_acp formula in
% 'AcylACPsynthesis' reaction. Ignore for phosphatidate synthesis reactions
% specific for each of the lipid groups for which we have data. 

% first, compute acyl_acp formula
rxnIdx = find(strcmp('AcylACPsynthesis',PanGEM.rxns));
metsIdx = find(PanGEM.S(:,rxnIdx)<0); % LHS mets
mets = PanGEM.mets(metsIdx);
[Ematrix, elements] = computeElementalMatrix(PanGEM,mets,'generalFormula',true)
metCoefs = full(PanGEM.S(metsIdx,rxnIdx));
RHS_elements = abs(sum(Ematrix .* repmat(metCoefs,1,size(Ematrix,2)),1));
PanGEM.metFormulas{find(strcmp('Acyl_acp',PanGEM.mets))} = strcat('C',mat2str(RHS_elements(1)),'H',mat2str(RHS_elements(3)),'O',mat2str(RHS_elements(2)),'S',mat2str(RHS_elements(4)),'R');
PanGEM.metCharges(find(strcmp('Acyl_acp',PanGEM.mets))) = nansum(abs(metCoefs) .* PanGEM.metCharges(metsIdx));
[massImbalance, imBalancedMass, imBalancedCharge, imBalancedRxnBool, Elements, missingFormulaeBool, balancedMetBool] = checkMassChargeBalance(PanGEM);
imBalancedMass(rxnIdx)

% Now find all instances of 'acyl' in mets
substr = 'acyl';
acylIdx1 = find(~cellfun(@isempty,strfind(PanGEM.mets,substr)));
substr = 'Acyl';
acylIdx2 = find(~cellfun(@isempty,strfind(PanGEM.mets,substr)))
acylIdx = [acylIdx1;acylIdx2];


%% Compute formula and MW of each crude fraction

% DNA 
metToEval = 'DNA';
synthRxn = 'DNASynthesis'
metIdx = find(strcmp(metToEval,PanGEM.mets));
synthRxnIdx = find(strcmp(synthRxn,PanGEM.rxns));

constituentMets = [{'dATP'},{'dGTP'},{'dTTP'},{'dCTP'}];
for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},PanGEM.mets));
    BOFcoef(a) = full(PanGEM.S(constituentMetsIdx(a),synthRxnIdx));
end
constituentMets = [{'dAMP'},{'dGMP'},{'dTMP'},{'dCMP'}]; % a diphosphate is removed from each precursor
[Ematrix, elements] = computeElementalMatrix(PanGEM,constituentMets)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 30.973762 0];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr = 1000/metToEval_MW_total; % correction to 1000 g mol-1

tempModel = PanGEM;
tempModel.S(constituentMetsIdx,synthRxnIdx) = MWcorr * PanGEM.S(constituentMetsIdx,synthRxnIdx);

constituentMets = [{'dATP'},{'dGTP'},{'dTTP'},{'dCTP'}];
for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},tempModel.mets));
    BOFcoef(a) = full(tempModel.S(constituentMetsIdx(a),synthRxnIdx));
end
constituentMets = [{'dAMP'},{'dGMP'},{'dTMP'},{'dCMP'}]; % a diphosphate is removed from each precursor
[Ematrix, elements] = computeElementalMatrix(tempModel,constituentMets)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 30.973762 0];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

BOFcoef
sum(BOFcoef)
elements
metToEval_formula
PanGEM.metFormulas{metIdx} = strcat('C',mat2str(metToEval_formula(1)), ...
    'H',mat2str(metToEval_formula(4)),'O',mat2str(metToEval_formula(3)), ...
    'N',mat2str(metToEval_formula(2)),'P',mat2str(metToEval_formula(5)));
PanGEM.S(constituentMetsIdx,synthRxnIdx) = BOFcoef;
clear constituentMets constituentMetsIdx BOFcoef Ematrix elements stoichMat metToEval_formula metToEval_MW metToEval_MW_elements metToEval_MW_total metToEval MWcorr 


% RNA 
metToEval = 'RNA';
synthRxn = 'RNASynthesis'
metIdx = find(strcmp(metToEval,PanGEM.mets));
synthRxnIdx = find(strcmp(synthRxn,PanGEM.rxns));

constituentMets = [{'ATP'},{'GTP'},{'UTP'},{'CTP'}];
for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},PanGEM.mets));
    BOFcoef(a) = full(PanGEM.S(constituentMetsIdx(a),synthRxnIdx));
end
constituentMets = [{'AMP'},{'GMP'},{'UMP'},{'CMP'}]; % a diphosphate is removed from each precursor
[Ematrix, elements] = computeElementalMatrix(PanGEM,constituentMets)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 30.973762 0];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr = 1000/metToEval_MW_total; % correction to 1000 g mol-1

tempModel = PanGEM;
tempModel.S(constituentMetsIdx,synthRxnIdx) = MWcorr * PanGEM.S(constituentMetsIdx,synthRxnIdx);

constituentMets = [{'ATP'},{'GTP'},{'UTP'},{'CTP'}];
for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},tempModel.mets));
    BOFcoef(a) = full(tempModel.S(constituentMetsIdx(a),synthRxnIdx));
end
constituentMets = [{'AMP'},{'GMP'},{'UMP'},{'CMP'}]; % a diphosphate is removed from each precursor
[Ematrix, elements] = computeElementalMatrix(tempModel,constituentMets)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 30.973762 0];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

BOFcoef
sum(BOFcoef)
elements
metToEval_formula
PanGEM.metFormulas{metIdx} = strcat('C',mat2str(metToEval_formula(1)), ...
    'H',mat2str(metToEval_formula(4)),'O',mat2str(metToEval_formula(3)), ...
    'N',mat2str(metToEval_formula(2)),'P',mat2str(metToEval_formula(5)));
PanGEM.S(constituentMetsIdx,synthRxnIdx) = BOFcoef;
clear constituentMets constituentMetsIdx BOFcoef Ematrix elements stoichMat metToEval_formula metToEval_MW metToEval_MW_elements metToEval_MW_total metToEval MWcorr 

% Lipids 
metToEval = 'Lipid';
synthRxn = 'LipidSynthesis'
metIdx = find(strcmp(metToEval,PanGEM.mets));
synthRxnIdx = find(strcmp(synthRxn,PanGEM.rxns));

constituentMets = [{'Phosphatidylglycerol'},{'1_2_Diacyl_3_beta_D_galactosyl_sn_glycerol_MGDG'},{'Digalactosyl_diacylglycerol'},{'Sulfoquinovosyldiacylglycerol'}];
for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},PanGEM.mets));
    BOFcoef(a) = full(PanGEM.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(PanGEM,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 15.999 1.00784 30.973762 0 32.065];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr = 1000/metToEval_MW_total; % correction to 1000 g mol-1

tempModel = PanGEM;
tempModel.S(constituentMetsIdx,synthRxnIdx) = MWcorr * PanGEM.S(constituentMetsIdx,synthRxnIdx);

constituentMets = [{'Phosphatidylglycerol'},{'1_2_Diacyl_3_beta_D_galactosyl_sn_glycerol_MGDG'},{'Digalactosyl_diacylglycerol'},{'Sulfoquinovosyldiacylglycerol'}];
for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},tempModel.mets));
    BOFcoef(a) = full(tempModel.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(tempModel,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 15.999 1.00784 30.973762 0 32.065];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

BOFcoef
sum(BOFcoef)
elements
metToEval_formula

PanGEM.metFormulas{metIdx} = strcat('C',mat2str(metToEval_formula(1)), ...
    'H',mat2str(metToEval_formula(3)),'O',mat2str(metToEval_formula(2)), ...
    'P',mat2str(metToEval_formula(4)),'S',mat2str(metToEval_formula(6)));
PanGEM.S(constituentMetsIdx,synthRxnIdx) = BOFcoef;
clear constituentMets constituentMetsIdx BOFcoef Ematrix elements stoichMat metToEval_formula metToEval_MW metToEval_MW_elements metToEval_MW_total metToEval MWcorr 


% Protein 
metToEval = 'Protein';
synthRxn = 'ProteinSynthesis'
metIdx = find(strcmp(metToEval,PanGEM.mets));
synthRxnIdx = find(strcmp(synthRxn,PanGEM.rxns));

constituentMets = [{'Glycine'},{'L_Alanine'},{'L_Arginine'},{'L_Asparagine'}, ...
    {'L_Aspartate'},{'L_Cystine'},{'L_Glutamate'},{'L_Glutamine'}, ...
    {'L_Histidine'},{'L_Isoleucine'},{'L_Leucine'},{'L_Lysine'}, ...
    {'L_Methionine'},{'L_Phenylalanine'},{'L_Proline'},{'L_Selenomethionine'}, ...
    {'L_Serine'},{'L_Threonine'},{'L_Tryptophan'},{'L_Tyrosine'},{'L_Valine'}];

for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},PanGEM.mets));
    BOFcoef(a) = full(PanGEM.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(PanGEM,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 32.065 78.96];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr = 1000/metToEval_MW_total; % correction to 1000 g mol-1

tempModel = PanGEM;
tempModel.S(constituentMetsIdx,synthRxnIdx) = MWcorr * PanGEM.S(constituentMetsIdx,synthRxnIdx);

for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},tempModel.mets));
    BOFcoef(a) = full(tempModel.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(tempModel,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 32.065 78.96];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

BOFcoef
sum(BOFcoef)
elements
metToEval_formula
PanGEM.metFormulas{metIdx} = strcat('C',mat2str(metToEval_formula(1)), ...
    'H',mat2str(metToEval_formula(4)),'O',mat2str(metToEval_formula(3)), ...
    'N',mat2str(metToEval_formula(2)),'S',mat2str(metToEval_formula(5)), ...
    'Se',mat2str(metToEval_formula(6)));
PanGEM.S(constituentMetsIdx,synthRxnIdx) = BOFcoef;
clear constituentMets constituentMetsIdx BOFcoef Ematrix elements stoichMat metToEval_formula metToEval_MW metToEval_MW_elements metToEval_MW_total metToEval MWcorr 

% BioPool 
metToEval = 'BioPool';
synthRxn = 'BioPoolSynthesis'
metIdx = find(strcmp(metToEval,PanGEM.mets));
synthRxnIdx = find(strcmp(synthRxn,PanGEM.rxns));

constituentMets = [{'CDP_3_6_dideoxy_D_mannose'},{'Chorismate'},{'Glutathione'},{'Heme_O'}, ...
    {'S_Adenosyl_L_methionine'},{'Spermidine'}];


for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},PanGEM.mets));
    BOFcoef(a) = full(PanGEM.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(PanGEM,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 30.973762 32.065 55.845];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr = 1000/metToEval_MW_total; % correction to 1000 g mol-1

tempModel = PanGEM;
tempModel.S(constituentMetsIdx,synthRxnIdx) = MWcorr * PanGEM.S(constituentMetsIdx,synthRxnIdx);

for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},tempModel.mets));
    BOFcoef(a) = full(tempModel.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(tempModel,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 30.973762 32.065 55.845];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr
BOFcoef
sum(BOFcoef)
elements
metToEval_formula

PanGEM.metFormulas{metIdx} = strcat('C',mat2str(metToEval_formula(1)), ...
    'H',mat2str(metToEval_formula(4)),'O',mat2str(metToEval_formula(3)), ...
    'N',mat2str(metToEval_formula(2)),'P',mat2str(metToEval_formula(5)), ...
    'S',mat2str(metToEval_formula(6)),'Fe',mat2str(metToEval_formula(7)));
PanGEM.S(constituentMetsIdx,synthRxnIdx) = BOFcoef;
clear constituentMets constituentMetsIdx BOFcoef Ematrix elements stoichMat metToEval_formula metToEval_MW metToEval_MW_elements metToEval_MW_total metToEval MWcorr 

% CellWall
metToEval = 'CellWall';
synthRxn = 'CellWallSynthesis'
metIdx = find(strcmp(metToEval,PanGEM.mets));
synthRxnIdx = find(strcmp(synthRxn,PanGEM.rxns));

constituentMets = [{'di_trans_poly_cis_Undecaprenyl_diphosphate'},{'Lipid_A_disaccharide'},{'Crosslinked_peptideglycan'}];


for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},PanGEM.mets));
    BOFcoef(a) = full(PanGEM.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(PanGEM,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 30.973762];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr = 1000/metToEval_MW_total; % correction to 1000 g mol-1

tempModel = PanGEM;
tempModel.S(constituentMetsIdx,synthRxnIdx) = MWcorr * PanGEM.S(constituentMetsIdx,synthRxnIdx);

for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},tempModel.mets));
    BOFcoef(a) = full(tempModel.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(tempModel,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 30.973762];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr
BOFcoef
sum(BOFcoef)
elements
metToEval_formula

PanGEM.metFormulas{metIdx} = strcat('C',mat2str(metToEval_formula(1)), ...
    'H',mat2str(metToEval_formula(4)),'O',mat2str(metToEval_formula(3)), ...
    'N',mat2str(metToEval_formula(2)),'P',mat2str(metToEval_formula(5)));
PanGEM.S(constituentMetsIdx,synthRxnIdx) = BOFcoef;
clear constituentMets constituentMetsIdx BOFcoef Ematrix elements stoichMat metToEval_formula metToEval_MW metToEval_MW_elements metToEval_MW_total metToEval MWcorr 


% NB
metToEval = 'NB';
synthRxn = 'NBSynthesis'
metIdx = find(strcmp(metToEval,PanGEM.mets));
synthRxnIdx = find(strcmp(synthRxn,PanGEM.rxns));

constituentMets = [{'Adenine'},{'Adenosine'},{'Cytidine'},{'Cytosine'},...
    {'Deoxyadenosine'},{'Deoxycytidine'},{'Deoxyguanosine'},{'Deoxyribose'}, ...
    {'Deoxyuridine'},{'Guanine'},{'Guanosine'},{'Inosine'},{'Thymidine'}, ...
    {'Uracil'},{'Uridine'},{'Xanthosine'}];


for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},PanGEM.mets));
    BOFcoef(a) = full(PanGEM.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(PanGEM,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr = 1000/metToEval_MW_total; % correction to 1000 g mol-1

tempModel = PanGEM;
tempModel.S(constituentMetsIdx,synthRxnIdx) = MWcorr * PanGEM.S(constituentMetsIdx,synthRxnIdx);

for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},tempModel.mets));
    BOFcoef(a) = full(tempModel.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(tempModel,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr
BOFcoef
sum(BOFcoef)
elements
metToEval_formula

PanGEM.metFormulas{metIdx} = strcat('C',mat2str(metToEval_formula(1)), ...
    'H',mat2str(metToEval_formula(4)),'O',mat2str(metToEval_formula(3)), ...
    'N',mat2str(metToEval_formula(2)));
PanGEM.S(constituentMetsIdx,synthRxnIdx) = BOFcoef;
clear constituentMets constituentMetsIdx BOFcoef Ematrix elements stoichMat metToEval_formula metToEval_MW metToEval_MW_elements metToEval_MW_total metToEval MWcorr 


% Pigments
metToEval = 'Pigments';
synthRxn = 'PigmentSynthesis'
metIdx = find(strcmp(metToEval,PanGEM.mets));
synthRxnIdx = find(strcmp(synthRxn,PanGEM.rxns));

constituentMets = [{'3Z_Phycocyanobilin'},{'3Z_Phycoerythrobilin'}, ...
    {'7_8_Dihydro_beta_carotene'},{'all_trans_zeta_Carotene'},...
    {'beta_Carotene'},{'delta_Carotene'},{'epsilon_Carotene'},{'gamma_Carotene'}, ...
    {'Lycopene'},{'Divinylchlorophyll_a'},{'Divinylchlorophyll_b'}, ...
    {'alpha_Carotene'},{'Zeaxanthin'}];


for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},PanGEM.mets));
    BOFcoef(a) = full(PanGEM.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(PanGEM,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 24.305];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr = 1000/metToEval_MW_total; % correction to 1000 g mol-1

tempModel = PanGEM;
tempModel.S(constituentMetsIdx,synthRxnIdx) = MWcorr * PanGEM.S(constituentMetsIdx,synthRxnIdx);

for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},tempModel.mets));
    BOFcoef(a) = full(tempModel.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(tempModel,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 24.305];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr
BOFcoef
sum(BOFcoef)
elements
metToEval_formula

PanGEM.metFormulas{metIdx} = strcat('C',mat2str(metToEval_formula(1)), ...
    'H',mat2str(metToEval_formula(4)),'O',mat2str(metToEval_formula(3)), ...
    'N',mat2str(metToEval_formula(2)),'Mg',mat2str(metToEval_formula(5)));
PanGEM.S(constituentMetsIdx,synthRxnIdx) = BOFcoef;
clear constituentMets constituentMetsIdx BOFcoef Ematrix elements stoichMat metToEval_formula metToEval_MW metToEval_MW_elements metToEval_MW_total metToEval MWcorr 


% VitaCofactors
metToEval = 'VitaCofactors';
synthRxn = 'VitaCofactorsSynthesis'
metIdx = find(strcmp(metToEval,PanGEM.mets));
synthRxnIdx = find(strcmp(synthRxn,PanGEM.rxns));

constituentMets = [{'10_Formyltetrahydrofolate'},{'5_10_Methylenetetrahydrofolate'}, ...
    {'5_Methyltetrahydrofolate'},{'Acetyl_CoA'},...
    {'Biotin'},{'CoA'},{'Cobamide_coenzyme'},{'FAD'}, ...
    {'Malonyl_CoA'},{'NAD'},{'NADH'}, ...
    {'NADP'},{'NADPH'},{'Riboflavin'},{'Tetrahydrofolate'},{'Thiamin_diphosphate'}];


for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},PanGEM.mets));
    BOFcoef(a) = full(PanGEM.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(PanGEM,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 30.973762 32.065 58.933195];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr = 1000/metToEval_MW_total; % correction to 1000 g mol-1

tempModel = PanGEM;
tempModel.S(constituentMetsIdx,synthRxnIdx) = MWcorr * PanGEM.S(constituentMetsIdx,synthRxnIdx);

for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},tempModel.mets));
    BOFcoef(a) = full(tempModel.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(tempModel,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 30.973762 32.065 58.933195];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr
BOFcoef
sum(BOFcoef)
elements
metToEval_formula

PanGEM.metFormulas{metIdx} = strcat('C',mat2str(metToEval_formula(1)), ...
    'H',mat2str(metToEval_formula(4)),'O',mat2str(metToEval_formula(3)), ...
    'N',mat2str(metToEval_formula(2)),'P',mat2str(metToEval_formula(5)), ...
    'S',mat2str(metToEval_formula(6)),'Co',mat2str(metToEval_formula(7)));
PanGEM.S(constituentMetsIdx,synthRxnIdx) = BOFcoef;
clear constituentMets constituentMetsIdx BOFcoef Ematrix elements stoichMat metToEval_formula metToEval_MW metToEval_MW_elements metToEval_MW_total metToEval MWcorr 



% Ions
metToEval = 'Ions';
synthRxn = 'IonSynthesis'
metIdx = find(strcmp(metToEval,PanGEM.mets));
synthRxnIdx = find(strcmp(synthRxn,PanGEM.rxns));

constituentMets = [{'Cadmium'},{'Calcium_cation'}, ...
    {'Chloride_ion'},{'Cobalt_ion'},...
    {'Copper'},{'Fe2'},{'Magnesium_cation'},{'Guanylyl_molybdenum_cofactor'}, ...
    {'K'},{'Sodium_cation'},{'Strontium_cation'}, ...
    {'Zn2'}];


for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},PanGEM.mets));
    BOFcoef(a) = full(PanGEM.S(constituentMetsIdx(a),synthRxnIdx));
end
PanGEM.metFormulas{constituentMetsIdx(8)} = 'C20H24N10O15P2S2Mo';
[Ematrix, elements] = computeElementalMatrix(PanGEM,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 30.973762 112.411 40.078 35.453 58.933195 63.546 55.845 24.305 32.065 95.95 39.0983 22.989769 87.62 65.38];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr = 1000/metToEval_MW_total; % correction to 1000 g mol-1

tempModel = PanGEM;
tempModel.S(constituentMetsIdx,synthRxnIdx) = MWcorr * PanGEM.S(constituentMetsIdx,synthRxnIdx);

for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},tempModel.mets));
    BOFcoef(a) = full(tempModel.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(tempModel,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 30.973762 112.411 40.078 35.453 58.933195 63.546 55.845 24.305 32.065 95.95 39.0983 22.989769 87.62 65.38];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr
BOFcoef
sum(BOFcoef)
elements
metToEval_formula

PanGEM.metFormulas{metIdx} = strcat('C',mat2str(metToEval_formula(1)), ...
    'H',mat2str(metToEval_formula(4)),'O',mat2str(metToEval_formula(3)), ...
    'N',mat2str(metToEval_formula(2)),'P',mat2str(metToEval_formula(5)), ...
    'S',mat2str(metToEval_formula(13)),'Cd',mat2str(metToEval_formula(6)), ...
    'Ca',mat2str(metToEval_formula(7)),'Cl',mat2str(metToEval_formula(8)), ...
    'Co',mat2str(metToEval_formula(9)),'Cu',mat2str(metToEval_formula(10)), ...
    'Fe',mat2str(metToEval_formula(11)),'Mg',mat2str(metToEval_formula(12)), ...
    'Mo',mat2str(metToEval_formula(14)),'K',mat2str(metToEval_formula(15)), ...
    'Na',mat2str(metToEval_formula(16)),'Sr',mat2str(metToEval_formula(17)), ...
    'Zn',mat2str(metToEval_formula(18)));
PanGEM.S(constituentMetsIdx,synthRxnIdx) = BOFcoef;
clear constituentMets constituentMetsIdx BOFcoef Ematrix elements stoichMat metToEval_formula metToEval_MW metToEval_MW_elements metToEval_MW_total metToEval MWcorr 


% Carbohydrate
metToEval = 'Carbohydrate';
synthRxn = 'CarbohydrateSynthesis'
metIdx = find(strcmp(metToEval,PanGEM.mets));
synthRxnIdx = find(strcmp(synthRxn,PanGEM.rxns));

constituentMets = [{'Glycogen'}];


for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},PanGEM.mets));
    BOFcoef(a) = full(PanGEM.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(PanGEM,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 15.999 1.00784];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr = 1000/metToEval_MW_total; % correction to 1000 g mol-1

tempModel = PanGEM;
tempModel.S(constituentMetsIdx,synthRxnIdx) = MWcorr * PanGEM.S(constituentMetsIdx,synthRxnIdx);

for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},tempModel.mets));
    BOFcoef(a) = full(tempModel.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(tempModel,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 15.999 1.00784];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr
BOFcoef
sum(BOFcoef)
elements
metToEval_formula

PanGEM.metFormulas{metIdx} = strcat('C',mat2str(metToEval_formula(1)), ...
    'H',mat2str(metToEval_formula(3)),'O',mat2str(metToEval_formula(2)));
PanGEM.S(constituentMetsIdx,synthRxnIdx) = BOFcoef;

clear constituentMets constituentMetsIdx BOFcoef Ematrix elements stoichMat metToEval_formula metToEval_MW metToEval_MW_elements metToEval_MW_total metToEval MWcorr 

% and now biomass
metToEval = 'biomass';
synthRxn = 'BIOMASSCRUDE'
metIdx = find(strcmp(metToEval,PanGEM.mets));
synthRxnIdx = find(strcmp(synthRxn,PanGEM.rxns));

constituentMets = [{'DNA'},{'RNA'}, ...
    {'Carbohydrate'},{'Pigments'},...
    {'Protein'},{'CellWall'},{'VitaCofactors'},{'Ions'}, ...
    {'NB'},{'BioPool'},{'Divinylchlorophyll_a'},{'Divinylchlorophyll_b'}, ...
    {'alpha_Carotene'},{'Zeaxanthin'}];


for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},PanGEM.mets));
    BOFcoef(a) = full(PanGEM.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(PanGEM,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 30.973762 24.305 32.065 78.96 58.933195 112.411 40.078 35.453 63.546 55.845 95.95 39.0983 22.989769 87.62 65.38];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr = 1000/metToEval_MW_total; % correction to 1000 g mol-1

tempModel = PanGEM;
tempModel.S(constituentMetsIdx,synthRxnIdx) = MWcorr * PanGEM.S(constituentMetsIdx,synthRxnIdx);

for a = 1:numel(constituentMets)
    constituentMetsIdx(a) = find(strcmp(constituentMets{a},tempModel.mets));
    BOFcoef(a) = full(tempModel.S(constituentMetsIdx(a),synthRxnIdx));
end
[Ematrix, elements] = computeElementalMatrix(tempModel,constituentMets,'generalFormula',true)

stoichMat = Ematrix .* repmat(-BOFcoef',1,numel(elements));
metToEval_formula = sum(stoichMat,1);
elementMW = [12.011 14.0067 15.999 1.00784 30.973762 24.305 32.065 78.96 58.933195 112.411 40.078 35.453 63.546 55.845 95.95 39.0983 22.989769 87.62 65.38];
metToEval_MW_elements = metToEval_formula .* elementMW;
metToEval_MW_total = sum(metToEval_MW_elements);

MWcorr
BOFcoef
sum(BOFcoef)
elements
metToEval_formula

PanGEM.metFormulas{metIdx} = strcat('C',mat2str(metToEval_formula(1)), ...
    'H',mat2str(metToEval_formula(4)),'O',mat2str(metToEval_formula(3)), ...
    'N',mat2str(metToEval_formula(2)),'P',mat2str(metToEval_formula(5)), ...
    'S',mat2str(metToEval_formula(7)),'Se',mat2str(metToEval_formula(8)), ...
    'Cd',mat2str(metToEval_formula(10)), 'Ca',mat2str(metToEval_formula(11)),'Cl',mat2str(metToEval_formula(12)), ...
    'Co',mat2str(metToEval_formula(9)),'Cu',mat2str(metToEval_formula(13)), ...
    'Fe',mat2str(metToEval_formula(14)),'Mg',mat2str(metToEval_formula(6)), ...
    'Mo',mat2str(metToEval_formula(15)),'K',mat2str(metToEval_formula(16)), ...
    'Na',mat2str(metToEval_formula(17)),'Sr',mat2str(metToEval_formula(18)), ...
    'Zn',mat2str(metToEval_formula(19)));
PanGEM.S(constituentMetsIdx,synthRxnIdx) = BOFcoef;
clear constituentMets constituentMetsIdx BOFcoef Ematrix elements stoichMat metToEval_formula metToEval_MW metToEval_MW_elements metToEval_MW_total metToEval MWcorr 


sol = optimizeCbModel(PanGEM,[],'one');
col1 = PanGEM.rxns(find(strcmp('Exchange',PanGEM.subSystems)));
col2 = sol.x(find(strcmp('Exchange',PanGEM.subSystems)));
table(col1,col2)


%% Add rev field

for a = 1:numel(PanGEM.rxns)
    isRev(a) = PanGEM.lb(a) < 0 & PanGEM.ub(a) > 0;
end
PanGEM.rev = isRev';

%% Update PanGEM description
PanGEM.description = 'Prochlorococcus pangenome PanGEM - annotated and processed';

%% Save annotated PanGEM
PanGEM.description = 'Step 5 PanGEM'

fileName = 'data/models/Pro_PanGEM.mat';
save(fileName,'PanGEM');



