function [near_fieldT, near_fieldQ] = cal_near_field(mesh,params_FMM,gq_number,voxel_size)
%% Surface parameter
fprintf('Near field calculation begins.\n');
OT = params_FMM.OcTree;
Mask = OT.Mask;
Neighbor_list = params_FMM.Neighbor_list;
matrix_size = size(Mask);
TR = triangulation(mesh.faces,double(mesh.vertices));
nv = faceNormal(TR);
nf = size(mesh.faces,1);
[Xq, W] = quadrature(mesh, gq_number);

[Y,X,Z]=meshgrid(-(matrix_size(2)-1)/2:(matrix_size(2)-1)/2,...
    -(matrix_size(1)-1)/2:(matrix_size(1)-1)/2,...
    -(matrix_size(3)-1)/2:(matrix_size(3)-1)/2);

X = X*voxel_size(1);
Y = Y*voxel_size(2);
Z = Z*voxel_size(3);

leaf_depth = max(OT.BinDepths);
cell_leaf = find(OT.BinDepths==leaf_depth);

%% Precompute the number of near-field integrations
max_nonzeros = 0;
for i = 1:length(cell_leaf)
    leaf_i = cell_leaf(i);
    Index_Range = OT.BinIndex_Range(leaf_i,:);
    cell_mask = Mask(Index_Range(1):Index_Range(4), ...
        Index_Range(2):Index_Range(5), ...
        Index_Range(3):Index_Range(6));
    nnz_cell_mask = nnz(cell_mask);
    if nnz_cell_mask == 0
        continue;
    end
    neighbor_cell = [Neighbor_list(leaf_i,:) leaf_i];
    neighbor_cell = neighbor_cell(~isnan(neighbor_cell));
    neighbor_point = ismember(OT.PointBins, neighbor_cell);
    nnz_neighbor_point = nnz(neighbor_point);
    if nnz_neighbor_point == 0
        continue;
    end
    max_nonzeros = max_nonzeros + nnz_cell_mask * nnz_neighbor_point;
end

%% Initialize arrays
row_indices = zeros(max_nonzeros, 1);
col_indices = zeros(max_nonzeros, 1);
values_T = zeros(max_nonzeros, 1);
values_Q = zeros(max_nonzeros, 1);
current_index = 1;

Xq2 = reshape(Xq,[gq_number,size(mesh.faces,1),3]);
W2 = reshape(W,[gq_number,size(mesh.faces,1)]);
%% QUADRATURE
total_iterations = length(cell_leaf);
fprintf('QUADRATURE progress: %.0f%% completed, elapsed time: %d s\n', 0, 0);
tic
for leaf_i = 1:length(cell_leaf)
    % Check and display progress
    if mod(leaf_i, round(total_iterations * 0.5)) == 0 % Update every 50%
        elapsed_time = toc; % Get elapsed time in seconds
        fprintf('QUADRATURE progress: %.0f%% completed, elapsed time: %d s\n' ...
            , leaf_i/total_iterations*100, round(elapsed_time));
    end
    Index_Range = OT.BinIndex_Range(cell_leaf(leaf_i),:);
    cell_mask = Mask(Index_Range(1):Index_Range(4),Index_Range(2):Index_Range(5),Index_Range(3):Index_Range(6));
    nnz_cell_mask = nnz(cell_mask);
    if nnz_cell_mask==0
        continue
    end
    neighbor_cell = [Neighbor_list(cell_leaf(leaf_i),:) cell_leaf(leaf_i)];
    neighbor_cell = neighbor_cell(~isnan(neighbor_cell));
    neighbor_point = ismember(OT.PointBins,neighbor_cell);
    nnz_neighbor_point = nnz(neighbor_point);
    if nnz_neighbor_point==0
        continue
    end
    % Find the linear indices of values equal to 1 in cell_mask
    linear_indices = find(cell_mask);

    % Convert linear indices to subscripts (local coordinates)
    [cell_x, cell_y, cell_z] = ind2sub(size(cell_mask), linear_indices);

    % Convert local coordinates to global coordinates (a, b, c)
    a = Index_Range(1) - 1 + cell_x;
    b = Index_Range(2) - 1 + cell_y;
    c = Index_Range(3) - 1 + cell_z;

    % Calculate the linear indices in the entire mask
    indices_in_mask = sub2ind(size(Mask), a, b, c);
    u_coor_phy = [X(indices_in_mask), Y(indices_in_mask), Z(indices_in_mask)];

    % Transform dimensions
    u_coor_phy = repmat(u_coor_phy,[1,1,nnz_neighbor_point*gq_number]);
    u_coor_phy = permute(u_coor_phy,[1,3,2]);
    nv_i = repmat(nv(neighbor_point,:),1,1,gq_number);
    nv_i = permute(nv_i,[3,1,2]);
    nv_i = reshape(nv_i,[nnz_neighbor_point*gq_number,3]);
    Xq_i = Xq2(:,neighbor_point,:);
    W_i = W2(:,neighbor_point);
    Xq_i = reshape(Xq_i,[nnz_neighbor_point*gq_number,3]);
    W_i = reshape(W_i,[nnz_neighbor_point*gq_number,1]);

    Xq_i = repmat(Xq_i,[1,1,nnz_cell_mask]);
    Xq_i = permute(Xq_i,[3,1,2]);
    W_i = repmat(W_i,[1,nnz_cell_mask]);
    W_i = permute(W_i,[2,1]);
    nv_i = repmat(nv_i,[1,1,nnz_cell_mask]);
    nv_i = permute(nv_i,[3,1,2]);
    r_square = @(x) sum((x-u_coor_phy).^2,3);

    nf_T = gather(1./(4*pi)*sum(reshape(W_i.*r_square(Xq_i).^-0.5,[nnz_cell_mask,gq_number,nnz_neighbor_point]),2));
    nf_Q = gather(1/(-4*pi).*sum(reshape(W_i.*(r_square(Xq_i).^(-3/2) ...
        .*sum((Xq_i-u_coor_phy).*nv_i,3)),[nnz_cell_mask,gq_number,nnz_neighbor_point]),2));

    num_voxels = length(indices_in_mask);
    neighbor_point_indices = find(neighbor_point);
    num_neighbors = length(neighbor_point_indices);

    nf_T = reshape(nf_T, [num_voxels, num_neighbors]);
    nf_Q = reshape(nf_Q, [num_voxels, num_neighbors]);

    % Generate row and column indices
    [rows, cols] = ndgrid(indices_in_mask, neighbor_point_indices);
    num_entries = numel(rows);

    idx_range = current_index : current_index + num_entries - 1;

    row_indices(idx_range) = rows(:);
    col_indices(idx_range) = cols(:);
    values_T(idx_range) = nf_T(:);
    values_Q(idx_range) = nf_Q(:);

    current_index = current_index + num_entries;
end

% Trim arrays
row_indices = row_indices(1:current_index - 1);
col_indices = col_indices(1:current_index - 1);
values_T = values_T(1:current_index - 1);
values_Q = values_Q(1:current_index - 1);

% Create sparse matrix index mapping
max_row = numel(Mask);
max_col = nf;
sparse_idx = sparse(double(row_indices(1:current_index - 1)), double(col_indices(1:current_index - 1)), ...
                    1:(current_index - 1), max_row, max_col);

%% Analytical computation of near-singular integrals
fprintf('Analytical intergral progress: %.0f%% completed, elapsed time: %d s\n', 0, 0);
point_indices_per_triangle = nearpoint_prer_triangle(params_FMM,mesh,voxel_size,8);
tic
for t = 1:nf
    % Check and display progress
    if mod(t, round(nf * 0.25)) == 0 % Update every 25%
        elapsed_time = toc; % Get elapsed time in seconds
        fprintf('Analytical intergral progress: %.0f%% completed, elapsed time: %d s\n' ...
            , t/round(nf)*100, round(elapsed_time));
    end
    % Get the indices of P_i associated with the t-th triangle
    idx_in_triangle = point_indices_per_triangle{t};

    % Skip if there are no points
    if isempty(idx_in_triangle)
        continue;
    end   
    % Get the vertex indices of the triangle
    verts_idx = mesh.faces(t, :);
    % Get the coordinates of the vertices
    V1 = mesh.vertices(verts_idx(1), :);
    V2 = mesh.vertices(verts_idx(2), :);
    V3 = mesh.vertices(verts_idx(3), :);
    % Normal vector of the current triangle
    nv_t = nv(t, :);
    %near-singular integrals
    [T, Q] = compute_near_singular_integral(idx_in_triangle,X,Y,Z,V1,V2,V3,nv_t);

    row_indices_t = idx_in_triangle;
    col_indices_t = repmat(t, numel(idx_in_triangle), 1);

    % Find index positions from the sparse matrix
    idx_ranges = full(sparse_idx(sub2ind([max_row, max_col], row_indices_t, col_indices_t)));

    % Find valid indices
    valid_indices = idx_ranges > 0;
    idx_ranges_valid = idx_ranges(valid_indices);
    T_valid = T(valid_indices);
    Q_valid = Q(valid_indices);

    % Replace values in values_T and values_Q
    values_T(idx_ranges_valid) = T_valid;
    values_Q(idx_ranges_valid) = Q_valid;
end

%% Create sparse matrices
near_fieldT = sparse(row_indices, col_indices, values_T, numel(Mask), nf);
near_fieldQ = sparse(row_indices, col_indices, values_Q, numel(Mask), nf);
fprintf('Near field calculation completed.\n');
end