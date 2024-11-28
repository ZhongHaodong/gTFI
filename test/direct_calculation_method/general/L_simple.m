function [L, mask] = L_simple(mesh, gq_number, voxel_size, mask)
%   [L, mask] = L_simple(mesh, gq_number, voxel_size, mask)
%   output:
%   L - the matrix form of L
%   mask - new mask
%
%   input:
%   mesh - triangular surface mesh
%          - mesh.vertices: Vertex coordinates (N x 3 matrix)
%          - mesh.faces: Triangular mesh faces (M x 3 matrix)
%   gq_number - the number of Gaussian quadrature nodes
%   voxel_size - voxel size in physical units
%   mask - binary mask for the domain
%
%   Created by Haodong zhong (zhonghaodonghy@outlook.com) on 2022.06
%   Last modified by Haodong zhong on 2024.11.25

%% Surface parameter
matrix_size = size(mask);

TR = triangulation(mesh.faces, mesh.vertices);
nf = size(mesh.faces, 1); % Number of faces
nv = faceNormal(TR); % Unit outward normal vector

%% Calculate the matrix D2N used to convert Neumann boundary conditions to Dirichlet boundary conditions
D2N = computeD2N(mesh, 64);

%% Precondition
[idx_x, idx_y, idx_z] = ind2sub(size(mask), find(mask));
u_coor = [idx_x, idx_y, idx_z]; % Extract the coordinates of non-zero points in the mask

[Xq, W] = quadrature(mesh, gq_number);
Xq = gpuArray(Xq);
W = gpuArray(W);
nv2 = repmat(nv, 1, 1, gq_number);
nv2 = permute(nv2, [3, 1, 2]);
nv2 = reshape(nv2, [size(mesh.faces, 1) * gq_number, 3]);

%% Calculate T1
[Y, X, Z] = meshgrid(-(matrix_size(2)-1)/2:(matrix_size(2)-1)/2, ...
                     -(matrix_size(1)-1)/2:(matrix_size(1)-1)/2, ...
                     -(matrix_size(3)-1)/2:(matrix_size(3)-1)/2);

X = X * voxel_size(1);
Y = Y * voxel_size(2);
Z = Z * voxel_size(3);

T1 = zeros([size(u_coor, 1), nf]);

total_iterations = size(u_coor, 1);
fprintf('T` progress: %.0f%% completed, elapsed time: %d minutes\n', 0, 0);
tic
for i = 1:total_iterations
    % Check and output progress
    elapsed_time = toc; % Elapsed time in seconds
    elapsed_minutes = floor(elapsed_time / 60); % Minutes
    elapsed_seconds = floor(mod(elapsed_time, 60)); % Seconds
    if mod(i, round(total_iterations * 0.2)) == 0 % Update every 20%
        fprintf('T` progress: %.0f%% completed, elapsed time: %d:%02d\n', i/total_iterations*100, elapsed_minutes, elapsed_seconds);
    end
    a = u_coor(i, 1);
    b = u_coor(i, 2);
    c = u_coor(i, 3);
    u_coor_phy = gpuArray([X(a, b, c), Y(a, b, c), Z(a, b, c)]);
    r_square = @(x) sum((x - u_coor_phy).^2, 2);

    fT = @(x) 1 ./ ((4 * pi) * r_square(x).^0.5);
    I0 = W .* fT(Xq);
    I1 = reshape(I0, [gq_number, size(mesh.faces, 1)]);
    I2 = permute(I1, [2, 1]);
    T1(i, :) = sum(I2, 2);
end

L = T1 * D2N;
clear D2N T1

%% Calculate Q1
Q1 = zeros([size(u_coor, 1), nf]);

fprintf('Q` progress: %.0f%% completed, elapsed time: %d minutes\n', 0, 0);
tic
for i = 1:total_iterations
    % Check and output progress
    elapsed_time = toc; % Elapsed time in seconds
    elapsed_minutes = floor(elapsed_time / 60); % Minutes
    elapsed_seconds = floor(mod(elapsed_time, 60)); % Seconds
    if mod(i, round(total_iterations * 0.2)) == 0 % Update every 20%
        fprintf('Q` progress: %.0f%% completed, elapsed time: %d:%02d\n', i/total_iterations*100, elapsed_minutes, elapsed_seconds);
    end
    a = u_coor(i, 1);
    b = u_coor(i, 2);
    c = u_coor(i, 3);
    u_coor_phy = gpuArray([X(a, b, c), Y(a, b, c), Z(a, b, c)]);
    r_square = @(x) sum((x - u_coor_phy).^2, 2);

    fQ = @(x) 1 / (-4 * pi) .* r_square(x).^(-3/2) .* sum((x - u_coor_phy) .* nv2, 2);
    I0 = W .* fQ(Xq);
    I1 = reshape(I0, [gq_number, size(mesh.faces, 1)]);
    I2 = permute(I1, [2, 1]);
    Q1(i, :) = sum(I2, 2);
end

clear Xq W nv2 I0 I1 I2
L = L - Q1;
end
