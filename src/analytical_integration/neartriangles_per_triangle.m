function triangle_indices_per_triangle = neartriangles_per_triangle(params_FMM, mesh, distance)
% neartriangles_per_triangle: Assigns triangles within a specified distance to each triangle.
%
% This function computes the indices of triangles whose centroids are within a specified
% distance from each triangle's centroid. It uses a precomputed octree structure (OcTree) and
% neighborhood relationships to efficiently find nearby triangles.
%
% Inputs:
%   params_FMM - A structure containing precomputed octree data and neighbor relationships.
%                Expected fields include:
%       - params_FMM.OcTree.Centroid: Centroids of the triangles.
%       - params_FMM.OcTree.PointBins: Bin indices for each triangle centroid.
%       - params_FMM.Neighbor_list: Neighbors for each leaf bin (excluding itself), a 2D matrix.
%   mesh       - A structure representing the 3D mesh. Expected fields include:
%       - mesh.faces: (Nx3) Array of triangle vertex indices.
%       - mesh.vertices: (Mx3) Array of vertex coordinates.
%   distance   - The threshold distance to consider when finding nearby triangles.
%
% Outputs:
%   triangle_indices_per_triangle - A cell array where each cell corresponds to a triangle.
%                                   Each cell contains a list of triangle indices whose centroids
%                                   are within the specified distance from the current triangle's centroid.
%
%   Created on 2024.11.30

% Extract octree and neighbor information
OcTree = params_FMM.OcTree;
PointBins = OcTree.PointBins; % Bin index for each triangle centroid
Neighbor_list = params_FMM.Neighbor_list; % Neighbor bins for each bin (excluding itself), 2D matrix

% Number of triangles
n_faces = size(mesh.faces, 1);

% Initialize output cell array
triangle_indices_per_triangle = cell(n_faces, 1);

% Use precomputed centroids
centroids = OcTree.Centroid;

% Build a mapping from bins to triangle indices
max_bin_index = max(PointBins);
bin_to_triangle_indices = cell(max_bin_index, 1);

for t = 1:n_faces
    bin_idx = PointBins(t);
    if isempty(bin_to_triangle_indices{bin_idx})
        bin_to_triangle_indices{bin_idx} = t;
    else
        bin_to_triangle_indices{bin_idx}(end+1) = t;
    end
end

% Loop over each triangle
for t = 1:n_faces
    % Get the bin index of the current triangle
    bin_t = PointBins(t);
    
    % Get neighbor bins for the current bin, include itself
    neighbor_bins = Neighbor_list(bin_t, :); % Get the row corresponding to bin_t
    neighbor_bins = neighbor_bins(~isnan(neighbor_bins)); % Remove NaNs
    neighbor_bins = [neighbor_bins, bin_t]; % Include current bin
    neighbor_bins = unique(neighbor_bins); % Remove duplicates if any
    
    % Collect triangle indices from neighbor bins
    neighbor_triangle_indices = [];
    for idx = 1:length(neighbor_bins)
        bin_idx = neighbor_bins(idx);
        if bin_idx > length(bin_to_triangle_indices) || bin_idx < 1
            continue; % Skip invalid bin indices
        end
        neighbor_triangles = bin_to_triangle_indices{bin_idx};
        if ~isempty(neighbor_triangles)
            neighbor_triangle_indices = [neighbor_triangle_indices; neighbor_triangles(:)];
        end
    end
    
    % Remove duplicates and the current triangle itself
    neighbor_triangle_indices = unique(neighbor_triangle_indices);
    neighbor_triangle_indices(neighbor_triangle_indices == t) = [];
    
    % If no neighboring triangles, continue to next triangle
    if isempty(neighbor_triangle_indices)
        continue;
    end
    
    % Compute distances from the current centroid to neighboring centroids
    centroids_neighbors = centroids(neighbor_triangle_indices, :);
    centroid_t = centroids(t, :);
    
    distances = sqrt(sum((centroids_neighbors - centroid_t).^2, 2));
    
    % Find triangles within the specified distance
    within_distance = distances < distance;
    
    if any(within_distance)
        triangle_indices_per_triangle{t} = neighbor_triangle_indices(within_distance);
    end
end
end
