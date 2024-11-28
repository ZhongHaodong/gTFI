function mesh = adjustMeshOrientation(mesh, mask, voxel_size)
% Adjust the orientation of the mesh so that all normals point outward,
% using the local geometric center of the mask around each face.
%
% Input:
%   mesh - A struct with fields:
%          - mesh.vertices: Nx3 array of vertex coordinates
%          - mesh.faces: Mx3 array of triangular face indices
%   mask - A 3D binary array representing the region of interest
%   voxel_size - A 1x3 vector specifying the size of each voxel along x, y, z
%
% Output:
%   mesh - The updated mesh struct with adjusted face orientations

% Validate input
if ~isfield(mesh, 'vertices') || ~isfield(mesh, 'faces')
    error('Input mesh must contain "vertices" and "faces" fields.');
end

mask = logical(mask);

% Validate voxel_size
if nargin < 3 || numel(voxel_size) ~= 3
    error('voxel_size must be provided as a 1x3 vector.');
end

% Compute the physical coordinates of the mask voxels
matrix_size = size(mask);

% Create coordinate grids centered around zero
[Y_grid, X_grid, Z_grid] = meshgrid(...
    -(matrix_size(2)-1)/2 : (matrix_size(2)-1)/2, ...
    -(matrix_size(1)-1)/2 : (matrix_size(1)-1)/2, ...
    -(matrix_size(3)-1)/2 : (matrix_size(3)-1)/2);

% Scale the grids by voxel_size
X_grid = X_grid * voxel_size(1);
Y_grid = Y_grid * voxel_size(2);
Z_grid = Z_grid * voxel_size(3);

% Precompute the coordinates of all mask voxels
mask_indices = find(mask);
mask_x = X_grid(mask_indices);
mask_y = Y_grid(mask_indices);
mask_z = Z_grid(mask_indices);

mask_voxel_coords = [mask_x, mask_y, mask_z];

% Extract vertices and faces
vertices = mesh.vertices;
faces = mesh.faces;

% Loop through each face to check and adjust orientation
for i = 1:size(faces, 1)
    % Current face vertex indices
    face = faces(i, :);

    % Get the coordinates of the three vertices
    v1 = vertices(face(1), :);
    v2 = vertices(face(2), :);
    v3 = vertices(face(3), :);

    % Compute the face normal (not normalized)
    normal = cross(v2 - v1, v3 - v1);

    % Compute the face centroid
    face_center = (v1 + v2 + v3) / 3;

    % Compute the radius as the distance from the centroid to any vertex
    radius = norm(face_center - v1);

    % Adjust the radius if it's less than 3 times voxel_size(1)
    min_radius = 2.5 * voxel_size(1);
    if radius < min_radius
        radius = min_radius;
    end

    % Find mask voxels within the sphere centered at face_center
    distances_sq = sum((mask_voxel_coords - face_center).^2, 2);
    within_sphere = distances_sq <= radius^2;

    % Get the local mask voxel coordinates
    local_mask_coords = mask_voxel_coords(within_sphere, :);

    % Check if there are any mask voxels within the sphere
    if ~isempty(local_mask_coords)
        % Compute the local mask center
        local_mask_center = mean(local_mask_coords, 1);

        % Compute vector from local mask center to face center
        vector_to_face = face_center - local_mask_center;

        % Check the orientation of the normal
        outward_check = dot(normal, vector_to_face);

        % If the normal points inward, reverse the vertex order
        if outward_check < 0
            faces(i, :) = face([1 3 2]); % Swap the second and third vertices
        end
    end
end

% Update the mesh with adjusted faces
mesh.faces = faces;
end
