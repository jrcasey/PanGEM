function AA_Distribution = getAADistribution(AASequence);
% get the frequency of occurances of each amino acid from a sequence

% Inputs
% AASequence -      String. 

% Outputs
% AA_Distribution - Structure. A structure containing the absolute counts,
%                   frequencies, and names of each amino acid. 

% John R. Casey
% 20190510

AA_Codes = [{'A'}, {'R'}, {'N'}, {'D'}, {'C'}, {'E'}, {'Q'}, {'G'}, {'H'},...
    {'I'}, {'L'}, {'K'}, {'M'}, {'F'}, {'P'}, {'S'}, {'T'}, {'W'}, {'Y'}, {'V'}];
AA_Names = [{'L_Alanine'}, {'L_Arginine'}, {'L_Asparagine'}, {'L_Aspartate'},...
    {'L_Cysteine'}, {'L_Glutamate'}, {'L_Glutamine'}, {'Glycine'}, {'L_Histidine'},...
    {'L_Isoleucine'}, {'L_Leucine'}, {'L_Lysine'}, {'L_Methionine'}, ...
    {'L_Phenylalanine'}, {'L_Proline'}, {'L_Serine'}, {'L_Threonine'}, {'L_Tryptophan'},...
    {'L_Tyrosine'}, {'L_Valine'}]; % GEM mets


nSequence = numel(AASequence);
nAAs = numel(AA_Codes);

for i = 1 : nAAs
    tempAA = AA_Codes{i};
    tempHits = strfind(AASequence, tempAA);
    nHits(i) = numel(tempHits);
end

AA_Frequency = nHits./nSequence;

AA_Distribution.codes = AA_Codes;
AA_Distribution.names = AA_Names;
AA_Distribution.count = nHits;
AA_Distribution.frequency = AA_Frequency;

end

