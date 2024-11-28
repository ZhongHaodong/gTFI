function point_indices_per_triangle = nearpoint_prer_triangle(params_FMM, mesh,voxel_size)
% nearpoint_prer_triangle: Assigns points within a specified distance to each triangle.
%
% This function computes the indices of points near each triangle in a 3D
% box. Points are grouped by their proximity to the triangle surfaces,
% based on a precomputed octree structure (OT) and neighborhood relationships.
%
% Inputs:
%   params_FMM - A structure containing precomputed octree data and neighbor
%                relationships. Expected fields include:
%       - params_FMM.OcTree: Octree structure with bin information.
%       - params_FMM.Neighbor_list: Precomputed neighbors for each octree bin.
%   mesh       - A structure representing the 3D mesh. Expected fields include:
%       - mesh.elt: (Nx3) Array of triangle vertex indices.
%       - mesh.vtx: (Mx3) Array of vertex coordinates.
%
% Outputs:
%   point_indices_per_triangle - A cell array where each cell corresponds to a
%                                triangle. Each cell contains a list of point
%                                indices within a specified distance of the
%                                triangle.
%
%   Created by Haodong zhong (zhonghaodonghy@outlook.com) on 2024.06
%   Last modified by Haodong zhong on 2024.11.26

OT = params_FMM.OcTree;
Mask = OT.Mask;
Neighbor_list = params_FMM.Neighbor_list;
matrix_size = size(Mask);
[Y,X,Z]=meshgrid(-matrix_size(2)/2:(matrix_size(2)/2-1),...
    -matrix_size(1)/2:(matrix_size(1)/2-1),...
    -matrix_size(3)/2:(matrix_size(3)/2-1));

X = X*voxel_size(1);
Y = Y*voxel_size(2);
Z = Z*voxel_size(3);
distance = 8;
% Number of faces (triangles)
n_faces = size(mesh.faces, 1);
point_indices_per_triangle = cell(n_faces, 1);

% Precompute cell to point indices mapping
n_cells = size(OT.BinIndex_Range, 1);
cell_to_point_indices = cell(n_cells, 1);

for cell_idx = 1:n_cells
    Index_Range = OT.BinIndex_Range(cell_idx, :);
    cell_mask = Mask(Index_Range(1):Index_Range(4), ...
        Index_Range(2):Index_Range(5), ...
        Index_Range(3):Index_Range(6));
    if any(cell_mask(:))
        [x_local, y_local, z_local] = ind2sub(size(cell_mask), find(cell_mask));
        x_global = Index_Range(1) - 1 + x_local;
        y_global = Index_Range(2) - 1 + y_local;
        z_global = Index_Range(3) - 1 + z_local;
        idx_in_cell = sub2ind(size(Mask), x_global, y_global, z_global);
        cell_to_point_indices{cell_idx} = idx_in_cell;
    else
        cell_to_point_indices{cell_idx} = [];
    end
end

% Loop over each triangle
for t = 1:n_faces
    % Get neighbor cells for the current triangle
    cell_t = OT.PointBins(t);
    neighbor_cell = [Neighbor_list(cell_t, :) cell_t];
    neighbor_cell = neighbor_cell(~isnan(neighbor_cell));

    % Collect point indices from neighbor cells
    neighbor_point_indices = cell_to_point_indices(neighbor_cell);

    % Filter out empty cells
    non_empty_cells = ~cellfun(@isempty, neighbor_point_indices);
    neighbor_point_indices = neighbor_point_indices(non_empty_cells);

    % If no points, skip to next triangle
    if isempty(neighbor_point_indices)
        continue;
    end

    % Concatenate point indices
    idx_in_bbox = vertcat(neighbor_point_indices{:});
    idx_in_bbox = unique(idx_in_bbox);

    % Get coordinates of these points
    X_pts = X(idx_in_bbox);
    Y_pts = Y(idx_in_bbox);
    Z_pts = Z(idx_in_bbox);
    points_in_bbox = [X_pts, Y_pts, Z_pts];

    % Get triangle vertices
    verts_idx = mesh.faces(t, :);
    V0 = mesh.vertices(verts_idx(1), :);
    V1 = mesh.vertices(verts_idx(2), :);
    V2 = mesh.vertices(verts_idx(3), :);

    % Compute triangle bounding box with padding
    min_coords = min([V0; V1; V2], [], 1) - distance;
    max_coords = max([V0; V1; V2], [], 1) + distance;

    % Filter points within the bounding box
    within_bbox = all(points_in_bbox >= min_coords & points_in_bbox <= max_coords, 2);

    if ~any(within_bbox)
        continue;
    end

    % Points within bounding box
    points_in_bbox = points_in_bbox(within_bbox, :);
    idx_in_bbox = idx_in_bbox(within_bbox);

    % Compute distances from points to triangle
    distances = point_to_triangle_distance(points_in_bbox, V0, V1, V2);

    % Find points within the specified distance
    within_distance = distances < distance;

    if any(within_distance)
        point_indices_per_triangle{t} = idx_in_bbox(within_distance);
    end
end
end
