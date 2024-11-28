function distances = point_to_triangle_distance(P, V0, V1, V2)
% point_to_triangle_distance: Calculates the distance from points P to the triangle defined by V0, V1, V2.
%
% Inputs:
%   P  - An N x 3 array of points, where each row represents a point's coordinates.
%   V0 - A 1 x 3 vector representing the first vertex of the triangle.
%   V1 - A 1 x 3 vector representing the second vertex of the triangle.
%   V2 - A 1 x 3 vector representing the third vertex of the triangle.
%
% Outputs:
%   distances - An N x 1 array representing the shortest distance from each point in P to the triangle.

% Compute the edges of the triangle
E0 = V1 - V0; % Edge V0->V1
E1 = V2 - V0; % Edge V0->V2

% Compute the normal vector of the triangle's plane
Nvec = cross(E0, E1); % Normal vector
Nlen = norm(Nvec); % Length of the normal vector
if Nlen == 0
    % Degenerate triangle (zero normal vector length)
    distances = NaN(size(P, 1), 1);
    return;
end
Nvec = Nvec / Nlen; % Normalize the normal vector

% Vector from V0 to P
D = bsxfun(@minus, P, V0); % N x 3

% Compute the signed distance from points to the triangle's plane
dist_plane = D * Nvec'; % N x 1

% Project points onto the triangle's plane
P_proj = P - dist_plane * Nvec; % Projected point coordinates (N x 3)

% Compute the vector from V0 to the projected points
D_proj = bsxfun(@minus, P_proj, V0); % N x 3

% Compute dot products for the triangle's edges
E0_dot_E0 = dot(E0, E0);
E0_dot_E1 = dot(E0, E1);
E1_dot_E1 = dot(E1, E1);

E0_dot_D_proj = D_proj * E0'; % Projection onto E0
E1_dot_D_proj = D_proj * E1'; % Projection onto E1

% Compute barycentric coordinates (s, t)
denom = E0_dot_E0 * E1_dot_E1 - E0_dot_E1 * E0_dot_E1; % Denominator
s = (E0_dot_E1 * E1_dot_D_proj - E1_dot_E1 * E0_dot_D_proj) / denom; % Barycentric coordinate s
t = (E0_dot_E1 * E0_dot_D_proj - E0_dot_E0 * E1_dot_D_proj) / denom; % Barycentric coordinate t

% Determine whether points are inside the triangle
inside = (s >= 0) & (t >= 0) & ((s + t) <= 1);

% Initialize the distances array
distances = zeros(size(P, 1), 1);

% For points inside the triangle, the distance is the absolute value of the signed distance
distances(inside) = abs(dist_plane(inside));

% For points outside the triangle, compute the shortest distance to the edges or vertices
outside = ~inside; % Points outside the triangle
if any(outside)
    P_outside = P(outside, :); % Coordinates of points outside the triangle

    % Compute distances to the triangle's three edges
    dist_e0 = point_to_segment_distance(P_outside, V0, V1);
    dist_e1 = point_to_segment_distance(P_outside, V1, V2);
    dist_e2 = point_to_segment_distance(P_outside, V2, V0);

    % Compute distances to the triangle's three vertices
    dist_v0 = sqrt(sum(bsxfun(@minus, P_outside, V0).^2, 2));
    dist_v1 = sqrt(sum(bsxfun(@minus, P_outside, V1).^2, 2));
    dist_v2 = sqrt(sum(bsxfun(@minus, P_outside, V2).^2, 2));

    % Select the minimum distance
    distances_outside = min([dist_e0, dist_e1, dist_e2, dist_v0, dist_v1, dist_v2], [], 2);
    distances(outside) = distances_outside; % Update distances for points outside
end
end
