function [T, Q] = compute_triangle_centroid_interactions(triangle_indices, params_FMM, V1, V2, V3, nv_t)
% compute_triangle_centroid_interactions: Computes the interaction between an input triangle
% and the centroids of neighboring triangles.
%
% Inputs:
%   triangle_indices - Indices of neighboring triangles.
%   centroids        - Nx3 array of centroids for all triangles.
%   V1, V2, V3       - 1x3 vectors representing the vertices of the input triangle.
%   nv_t             - 1x3 normal vector of the input triangle.
%
% Outputs:
%   T                - Computed T values for the centroids of neighboring triangles.
%   Q                - Computed Q values for the centroids of neighboring triangles.
%
% Note:
%   This function computes the interaction between the input triangle and the centroids of
%   neighboring triangles, adapting the method used in compute_near_singular_integral.
%
% Created on 2024.11.30

% 计算输入三角形的边向量
E12 = V2 - V1; % 边 V1->V2
E23 = V3 - V2; % 边 V2->V3
E31 = V1 - V3; % 边 V3->V1

% 计算单位边向量
n12 = E12 / norm(E12);
n23 = E23 / norm(E23);
n31 = E31 / norm(E31);

% 获取邻近三角形的重心坐标
centroids = params_FMM.OcTree.Centroid;
P_i = centroids(triangle_indices, :); % [n_points, 3]
n_points = size(P_i, 1);

V1_mat = repmat(V1, n_points, 1);
V2_mat = repmat(V2, n_points, 1);
V3_mat = repmat(V3, n_points, 1);

r_V1_Pi = V1_mat - P_i; % Size [n_points, 3]
r_V2_Pi = V2_mat - P_i;
r_V3_Pi = V3_mat - P_i;

% Compute signs based on dot product
v = sum(r_V1_Pi .* nv_t, 2);
signs = sign(v); % Size [n_points, 1]

% Compute p_k = r_Vk_Pi × n_{k(k+1)}
p1 = cross(r_V1_Pi, repmat(n12, n_points, 1)); % Size [n_points, 3]
p2 = cross(r_V2_Pi, repmat(n23, n_points, 1));
p3 = cross(r_V3_Pi, repmat(n31, n_points, 1));

% Compute angles between p vectors
cos_theta_p1_p2 = sum(p1 .* p2, 2) ./ (vecnorm(p1, 2, 2) .* vecnorm(p2, 2, 2));
cos_theta_p1_p2 = max(min(cos_theta_p1_p2, 1), -1); % Clamp values
angle_p1_p2 = abs(acos(cos_theta_p1_p2)); % In radians

cos_theta_p2_p3 = sum(p2 .* p3, 2) ./ (vecnorm(p2, 2, 2) .* vecnorm(p3, 2, 2));
cos_theta_p2_p3 = max(min(cos_theta_p2_p3, 1), -1);
angle_p2_p3 = abs(acos(cos_theta_p2_p3));

cos_theta_p3_p1 = sum(p3 .* p1, 2) ./ (vecnorm(p3, 2, 2) .* vecnorm(p1, 2, 2));
cos_theta_p3_p1 = max(min(cos_theta_p3_p1, 1), -1);
angle_p3_p1 = abs(acos(cos_theta_p3_p1));

% Compute theta
theta = signs .* (2 * pi - (angle_p1_p2 + angle_p2_p3 + angle_p3_p1));
Q = -1 / (4 * pi) * theta;

%% Compute T
o1 = sum(p1 .* repmat(nv_t, n_points, 1), 2);
o2 = sum(p2 .* repmat(nv_t, n_points, 1), 2);
o3 = sum(p3 .* repmat(nv_t, n_points, 1), 2);

norm_r_V1_Pi = vecnorm(r_V1_Pi, 2, 2);
norm_r_V2_Pi = vecnorm(r_V2_Pi, 2, 2);
norm_r_V3_Pi = vecnorm(r_V3_Pi, 2, 2);

dot_r_V1_Pi_n12 = sum(r_V1_Pi .* repmat(n12, n_points, 1), 2);
dot_r_V2_Pi_n23 = sum(r_V2_Pi .* repmat(n23, n_points, 1), 2);
dot_r_V3_Pi_n31 = sum(r_V3_Pi .* repmat(n31, n_points, 1), 2);

dot_r_V2_Pi_n12 = sum(r_V2_Pi .* repmat(n12, n_points, 1), 2);
dot_r_V3_Pi_n23 = sum(r_V3_Pi .* repmat(n23, n_points, 1), 2);
dot_r_V1_Pi_n31 = sum(r_V1_Pi .* repmat(n31, n_points, 1), 2);

% Compute logarithmic terms
L1 = log((norm_r_V2_Pi+dot_r_V2_Pi_n12)./(norm_r_V1_Pi+dot_r_V1_Pi_n12));
L2 = log((norm_r_V3_Pi+dot_r_V3_Pi_n23)./(norm_r_V2_Pi+dot_r_V2_Pi_n23));
L3 = log((norm_r_V1_Pi+dot_r_V1_Pi_n31)./(norm_r_V3_Pi+dot_r_V3_Pi_n31));

% Handle NaN and Inf values
L1(isnan(L1) | isinf(L1)) = 0;
L2(isnan(L2) | isinf(L2)) = 0;
L3(isnan(L3) | isinf(L3)) = 0;

% Final computation of T
T = 1 / (4 * pi) * (o1 .* L1 + o2 .* L2 + o3 .* L3 - v .* theta);
end
