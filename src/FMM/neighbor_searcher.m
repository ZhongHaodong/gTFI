function [Neighbor_list, Interaction_list] = neighbor_searcher(params_FMM)
% neighbor_searcher: Identifies neighbors and interaction lists for octree bins.
%
% Inputs:
%   OT:Structure representing the OCtree, containing:
%   OT.BinCount: Total number of OCTree bins.
%   OT.BinDepths: Depth of each bin.
%   OT.BinIndex_Range: Index range of bins.
%   OT.BinParents: Parent index of each bin.
%
% Outputs:
%   Neighbor_list    - List of neighbors for each bin (up to 26 neighbors).
%   Interaction_list - List of interacting bins for each bin (up to 189 interactions).
%
%   Created by Haodong zhong (zhonghaodonghy@outlook.com) on 2023.06
%   Last modified by Haodong zhong on 2024.11.26
% Initialize output arrays
OT = params_FMM.OcTree;
leaf_depth = max(OT.BinDepths); % Maximum depth of the octree
Neighbor_list = NaN(OT.BinCount, 26, 'single'); % Preallocate for 26 neighbors
Interaction_list = NaN(OT.BinCount, 189, 'single'); % Preallocate for 189 interactions

% Compute neighbors for each bin
for cell_depth = leaf_depth:-1:1
    % Find all bins at the current depth level
    cell_ThisLevel = find(OT.BinDepths == cell_depth);

    for i = 1:length(cell_ThisLevel)
        % Current bin index
        current_bin = cell_ThisLevel(i);

        % Check overlap in x, y, z dimensions
        x_same = logical(sum(ismember(OT.BinIndex_Range(cell_ThisLevel, [1, 4]) + [-0.5, 0.5], ...
            OT.BinIndex_Range(current_bin, [1, 4]) + [-0.5, 0.5]), 2));
        y_same = logical(sum(ismember(OT.BinIndex_Range(cell_ThisLevel, [2, 5]) + [-0.5, 0.5], ...
            OT.BinIndex_Range(current_bin, [2, 5]) + [-0.5, 0.5]), 2));
        z_same = logical(sum(ismember(OT.BinIndex_Range(cell_ThisLevel, [3, 6]) + [-0.5, 0.5], ...
            OT.BinIndex_Range(current_bin, [3, 6]) + [-0.5, 0.5]), 2));

        % Identify neighbors by overlap in all three dimensions
        same = x_same .* y_same .* z_same;
        Neighbor = cell_ThisLevel(logical(same)); % Neighbor indices at this level
        Neighbor = Neighbor(~ismember(Neighbor, current_bin)); % Exclude self

        % Assign neighbors to the Neighbor_list
        Neighbor_list(current_bin, 1:length(Neighbor)) = Neighbor;
    end
end

% Compute interactions for each bin
for i = 1:OT.BinCount
    % Only compute interactions for bins with depth >= 2
    if OT.BinDepths(i) >= 2
        % Parent bin of the current bin
        parents = OT.BinParents(i);

        % Interaction bins: bins with parent in neighbor list of the parent
        Interaction0 = find(ismember(OT.BinParents, Neighbor_list(parents, :)));

        % Exclude direct neighbors from the interaction bins
        Interaction = Interaction0(~ismember(Interaction0, Neighbor_list(i, :)));

        % Assign interactions to the Interaction_list
        Interaction_list(i, 1:length(Interaction)) = Interaction;
    end
end
end
