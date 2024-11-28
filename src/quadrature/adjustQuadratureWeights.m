function W = adjustQuadratureWeights(mesh, weights)
    % adjustQuadratureWeights
    % Adjusts quadrature weights to all triangles in the mesh

    vertices = mesh.vertices;
    faces = mesh.faces;
    % Compute triangle areas (norm of cross product / 2)
    areas = area(vertices(faces(:, 1), :), vertices(faces(:, 2), :), vertices(faces(:, 3), :));

    % Reshape areas and weights for broadcasting
    numTriangles = size(areas, 1);
    numWeights = size(weights, 1);
    reshapedAreas = reshape(areas, 1, numTriangles, 1);
    reshapedWeights = reshape(weights, numWeights, 1, 1);
    
    % Adjust weights for all triangles
    adjustedWeights = reshape(reshapedAreas .* reshapedWeights, numTriangles * numWeights, 1);
    
    % Output the adjusted weights
    W = adjustedWeights;
end
