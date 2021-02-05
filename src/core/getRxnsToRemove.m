function [rxnsToRemove] = getRxnsToRemove(model,database, strain)

%% Get exclusion list for all strains
% This function retrieves indices in a panGEM model for those reactions
% which should be removed for reconstructing each strain.
% Inputs
% model       -       Structure. Pangenome GEM structure
% database    -       Structure. Database structure for a clade or group of sequences
% strain      -       String. The name of a strain (needs to match to the database)
% Outputs
% rxnsToRemove -      Vector of indices to remove


% John R. Casey 20171207
Dat = database;
PanGEM = model;

strainIdx = find(strcmp(strain,Dat.orgDatabase.StrainName));
targStrKO_Idx = find(Dat.PresenceAbsenceMatrix(:,strainIdx(1)));
targStrKO = Dat.uniqueKO(targStrKO_Idx);

for a = 1:numel(targStrKO) % there's a blank, so let's start indexing at 2
    tempInd = find(strcmp(targStrKO{a},PanGEM.genes)); % find the kth strain KO's which are present in PanGEM
    if isempty(tempInd)
        targStrKO_PanGEM{a} = []; % shouldn't be any empties but just in case
    else
        targStrKO_PanGEM{a} = tempInd;
    end
end

targStrKO_PanGEM2 = cell2mat(targStrKO_PanGEM); % format
tempRxnGeneMat = full(PanGEM.rxnGeneMat(:,targStrKO_PanGEM2)); % retrieve reactions associated with strain genes
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
