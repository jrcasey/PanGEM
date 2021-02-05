%% Step 3 - Reconstruction

% This step may take several weeks - months if done carefully. I refer the
% user to ref for a comprehensive guide to this process. Briefly, we begin
% with the Pangenome Assembly geneerated in Step 2, using the list of
% unique KO's annotated in the pangenome. From these, use KEGG's REST API
% to link KO's to their respective reactions and populate an excel model. I
% used the formatting accepted by RAVEN, and this format will be maintained
% in Step 4. For guidance, take a look at the excel model in data/models/

% There are some automated pipelines for generating GEMs but as far as I'm
% aware, there is still no substitute for manual curation! Womp womp. I
% would highly recommend following this protocol (below) and picking up a
% copy of Berhard Palsson's book on Constraint-Based Modeling... a must
% have!

% Heirendt, L., Arreckx, S., Pfau, T., Mendoza, S. N., Richelle, A.,
% Heinken, A., Haraldsdóttir, H. S., Wachowiak, J., Keating, S. M., Vlasov,
% V., Magnusdóttir, S., Ng, C. Y., Preciat, G., ?agare, A., Chan, S. H. J.,
% Aurich, M. K., Clancy, C. M., Modamio, J., Sauls, J. T., ? Fleming, R. M.
% T. (2019). Creation and analysis of biochemical constraint-based models
% using the COBRA Toolbox v.3.0. Nature Protocols, 14(3), 639?702.
% https://doi.org/10.1038/s41596-018-0098-2


% Once you've finished (congratulations!), save as an excel model for easy
% reference and Step 4 will import it for you.



