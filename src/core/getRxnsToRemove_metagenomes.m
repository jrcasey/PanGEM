function [rxnsToRemove] = getRxnsToRemove_metagenomes(model,database, sample, alpha)

%% Get exclusion list for all metagenome samples
% This function retrieves indices in a panGEM model for those reactions
% which should be removed for reconstructing each metagenome sample GEMs.
% Inputs
% model       -       Structure. Pangenome GEM structure
% database    -       Structure. Database structure for a clade or group of sequences
% sample      -       String. The name of metagenome sample
% alpha       -       String. The confidence level used to generate KO Matrix (alpha = 'alpha_90', 'alpha_95', or 'alpha_99'). Default is alpha_95.
% Outputs
% rxnsToRemove -      Vector of indices to remove


% John R. Casey 20171207
% Tatsuro Tanioka 20220725

Dat = database;
PanGEM = model;

sampleIdx = find(strcmp(sample,Dat.SampleID));
% targSampleKO_Idx = find(Dat.PresenceAbsenceMatrix(:,sampleIdx(1)));
% targSampleKO = Dat.uniqueKO(targSampleKO_Idx);

if alpha == "alpha_90";
    targSampleKO_Idx = find(Dat.PresenceAbsenceMatrix_alpha90(:,sampleIdx(1)));
elseif alpha == "alpha_95";
    targSampleKO_Idx = find(Dat.PresenceAbsenceMatrix_alpha95(:,sampleIdx(1)));
elseif alpha == "alpha_99";
    targSampleKO_Idx = find(Dat.PresenceAbsenceMatrix_alpha99(:,sampleIdx(1)));
else
    targSampleKO_Idx = find(Dat.PresenceAbsenceMatrix_alpha95(:,sampleIdx(1)));
end

targSampleKO = Dat.uniqueKO(targSampleKO_Idx);

for a = 1:numel(targSampleKO) % there's a blank, so let's start indexing at 2
    tempInd = find(strcmp(targSampleKO{a},PanGEM.genes)); % find the kth strain KO's which are present in PanGEM
    if isempty(tempInd)
        targSampleKO_PanGEM{a} = []; % shouldn't be any empties but just in case
    else
        targSampleKO_PanGEM{a} = tempInd;
    end
end

targSampleKO_PanGEM2 = cell2mat(targSampleKO_PanGEM); % format
tempRxnGeneMat = full(PanGEM.rxnGeneMat(:,targSampleKO_PanGEM2)); % retrieve reactions associated with strain genes
[rowI, colI] = find(tempRxnGeneMat); % find indices of non-zero entries
targStrRxnsInd = sort(unique(rowI)); % get indices of reactions
rxnsToRemove = setdiff(1:numel(PanGEM.rxns),targStrRxnsInd); % retrieve reactions in PanGEM but not in the kth strain

% get reactions in PanGEM which need to be retained
columnSums = sum(full(PanGEM.rxnGeneMat),2);
retainRxnInd = find(columnSums==0);

rxnsToRemoveNew = setdiff(rxnsToRemove,retainRxnInd);


rxnsToRemove2 = PanGEM.rxns(rxnsToRemoveNew); % get the reaction ID's


rxnsToRemove = rxnsToRemove2;


end
