%   Adjust triangular faces to ensure their normals are not parallel to the z-axis  
%   Created by Haodong Zhong on 2023.07  
%   Last modified by Haodong Zhong on 2024.11.25

function mesh = adjust_mesh(mesh)
TR = triangulation(mesh.faces, double(mesh.vertices));

% Compute face normals
nv = faceNormal(TR);

% Extract vertices and faces
vertices = mesh.vertices;
faces = mesh.faces;

% Define a small threshold to determine if a vector is parallel to the z-axis
threshold = 1e-3;
parallel_indices = find(abs(nv(:, 1)) < threshold & abs(nv(:, 2)) < threshold);
fprintf('Found %d face normals parallel to the z-axis.\n', length(parallel_indices));

% Create a random adjustment function
adjust_z = @(z) z + (rand(1) - 0.5) * 0.15;

% Adjust the faces with normals parallel to the z-axis
for i = 1:length(parallel_indices)
    idx = parallel_indices(i);
    face = faces(idx, :);

    % Fine-tune the z-coordinates of vertices until the normals are no longer parallel to the z-axis
    success = false;
    max_attempts = 1000; % Maximum number of attempts to avoid an infinite loop
    attempts = 0;

    while ~success && attempts < max_attempts
        attempts = attempts + 1;

        % Backup original vertex coordinates
        original_vertices = vertices(face, :);

        % Adjust the z-coordinate of one vertex from the face
        for j = 1:3
            vertices(face(j), 3) = adjust_z(vertices(face(j), 3));
        end

        % Recompute the normals of all related faces
        success = true;
        for j = 1:length(face)
            related_faces = find(any(faces == face(j), 2));
            for k = 1:length(related_faces)
                f_idx = related_faces(k);
                v1 = vertices(faces(f_idx, 1), :);
                v2 = vertices(faces(f_idx, 2), :);
                v3 = vertices(faces(f_idx, 3), :);
                normal = cross(v2 - v1, v3 - v1);
                normal = normal / norm(normal); % Normalize the vector

                % Update nv
                nv(f_idx, :) = normal;

                % Check if the normal is parallel to the z-axis
                if abs(normal(1)) < threshold && abs(normal(2)) < threshold
                    success = false;
                    break;
                end
            end
            if ~success
                break;
            end
        end

        % Restore the original vertex coordinates if the adjustment fails
        if ~success
            vertices(face, :) = original_vertices;
        end
    end

    if attempts == max_attempts
        warning('Maximum number of attempts reached; some normals may still be parallel to the z-axis.');
    end
end

% Update the mesh vertices and faces
mesh.vertices = vertices;
mesh.faces = faces;
end




