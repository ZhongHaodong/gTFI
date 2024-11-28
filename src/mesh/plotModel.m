function plotModel(mesh)
% plotModel - Plots the triangular mesh and face normals
%
% INPUT:
%   mesh - Triangular surface mesh
%          - mesh.vertices: Vertex coordinates (N x 3 matrix)
%          - mesh.faces: Triangular mesh faces (M x 3 matrix)

    % Calculate barycentric coordinates for each face
    barycentric = (mesh.vertices(mesh.faces(:,1),:) + ...
                   mesh.vertices(mesh.faces(:,2),:) + ...
                   mesh.vertices(mesh.faces(:,3),:)) ./ 3;

    % Compute face normals
    TR = triangulation(mesh.faces, double(mesh.vertices));
    nv = faceNormal(TR);

    % Plot the mesh
    patch('Faces', mesh.faces, 'Vertices', mesh.vertices, ...
          'FaceAlpha', 0.2, 'LineStyle', ':', 'EdgeColor', 'blue');
    daspect([1 1 1]);
    view(3);
    axis tight;
    set(gca, 'YTick', -1:1:1);
    camlight;
    hold on;

    % Plot face normals
    nv = nv*2;
    quiver3(barycentric(:,1), barycentric(:,2), barycentric(:,3), ...
            nv(:,1), nv(:,2), nv(:,3), 0.5, 'color', 'r');
end
