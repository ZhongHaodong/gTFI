function [X, W] = quadrature(mesh, gq_number)
% Compute quadrature points and weights for a given triangular mesh
%
% INPUT:
%   mesh      - Struct containing:
%                 - mesh.vertices: Vertex coordinates (N x 3 matrix)
%                 - mesh.faces: Triangular mesh faces (M x 3 matrix)
%   gq_number - Number of Gaussian quadrature points, 7 or 64
%
% OUTPUT:
%   X         - Quadrature points (Q x 3 matrix, where Q = M * gq_number)
%   W         - Quadrature weights (Q x 1 vector)

    % Reference quadrature
    [x, w] = quadrature_ut(gq_number);
    
    % Quadrature coordinates
    X = zeros(size(x,1) * size(mesh.faces, 1), size(mesh.vertices, 2));
    for j = 1:size(x,1)
        idx = (j:size(x,1):size(x,1)*size(mesh.faces,1))';
        X(idx, :) = (1 - x(j,1) - x(j,2)) * mesh.vertices(mesh.faces(:,1), :) ...
                  + x(j,1) * mesh.vertices(mesh.faces(:,2), :) ...
                  + x(j,2) * mesh.vertices(mesh.faces(:,3), :);
    end
    
    % Quadrature weights
    W = adjustQuadratureWeights(mesh, w);
end

