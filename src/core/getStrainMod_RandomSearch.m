function [reducedMod, nRxns] = getStrainMod_RandoSearch(model,rxnsToRemove,nPermutations,minGrowth)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Algorithm 2 - Minimize network size. 
% Reconstructs strain models by searching for the minimum number of
% reactions required to synthesize biomass. The algorithm initializes with
% a list of reactions that should be removed from PanGEM, then tries to
% remove them one by one. If a deletion is lethal then it is kept. The
% order matters, so a number of random permutations is attempted and tbe
% smallest network is selected. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% Dat               -   Assembly structure. (required)
% PanGEM            -   Pangenome scale model. (required)
% nPermutations     -   number of permutations. (opt, default = 100)
% minGrowth         -   Minimum biomass growth rate. (opt, default = 0.05)
% scaffoldCutoff    -   Number of scaffolds which define a 'high quality' 
%                       genome. (opt, default = 30)
% fileName          -   Directory and filename to store model structure.
%                       (opt, default is to not save)


% Outputs           
% targStrMod        -   Structure containing all minimal strain models
% nRxns2            -   Array (double) of network sizes for each strain for
%                       each permutation
% targStrGrowth     -   Array of growth rates for each strain
% targStrFluxes     -   Cell array of fluxes for each strain
% targStrShadow     -   Cell array of shadow prices for each strain

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version 1 - 20171030

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - John R. Casey
% jrcasey@hawaii.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for p = 1:nPermutations
        tempMod = model; % reset model to PanGEM for each permutation
        [rxnsToRemoveRand order(p,:)] = datasample(rxnsToRemove,numel(rxnsToRemove),'Replace',false); % generate a random permutation of reactionsToBeRemoved
        for i = 1:numel(rxnsToRemoveRand)
            tempMod2 = removeReactions(tempMod,rxnsToRemoveRand{i},true,true,false); % Sequentially remove each reaction
            tempSol = solveLP(tempMod2); % solve
            if abs(tempSol.f)>minGrowth % need to keep this above 0.04... there is some reaction which when removed drops the growth rate by a half. Let's go back and find out which one this is!
                tempMod = tempMod2; % if the solution is still viable, keep the change to the model. Otherwise, discard the change (next loop iteration will go back to the previous version without overwriting tempMod with tempMod2)
            end
        end
        PermuteMod(p) = tempMod; % all permutation models
        nRxns(p) = numel(tempMod.rxns); % number of reactions in each model.
        if abs(tempSol.f)> minGrowth % check for feasiblity again
            finalCheck(p) = 1; % assign true for feasible
        else finalCheck(p) = 0; % assign false for infeasible
        end
    end
    excludeInd = find(~finalCheck); % find the failed reconstructions
    nRxns(excludeInd) = 100000; % assign an arbitrarily large number to the corresponding nRxns
    [junk, MinI] = min(nRxns); % winning model index (feasible but with least number of reactions)
    reducedMod = PermuteMod(MinI); % winning model 
end


