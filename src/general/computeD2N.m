function D2N = computeD2N(mesh, gq_number)
%   D2N = computeD2N(mesh, gq_number)
%   output
%   D2N - the matrix used to 
%   convert Neumann boundary conditions to Dirichlet boundary conditions.
%
%   input
%   mesh - triangular surface mesh
%          - mesh.vertices: Vertex coordinates (N x 3 matrix)
%          - mesh.faces: Triangular mesh faces (M x 3 matrix)
%   gq_number - the number of Gaussian quadrature nodes
%
%   When using the code, please cite 
%
%   Created by Haodong zhong (zhonghaodonghy@outlook.com) on 2022.06
%   Last modified by Haodong zhong on 2024.11.25

%% Surface parameter
TR = triangulation(mesh.faces, double(mesh.vertices));
centroid = (mesh.vertices(mesh.faces(:,1),:) + ...
            mesh.vertices(mesh.faces(:,2),:) + ...
            mesh.vertices(mesh.faces(:,3),:)) ./ 3; % Triangle centroids
nf = size(mesh.faces, 1); % Number of faces
nv = faceNormal(TR); % Unit outward normal vector

%% Calculate D2N
tri_num = nf; % Number of triangles
[Xq, W] = quadrature(mesh, gq_number); % Quadrature points and weights
Xq = reshape(Xq, [gq_number, tri_num, 3]);
Xq = permute(Xq, [2, 3, 1]);
Xq = gpuArray(single(Xq)); % Transfer to GPU for faster computation
W = reshape(W, [gq_number, tri_num, 1]);
W = permute(W, [2, 3, 1]);
W = gpuArray(single(W)); % Quadrature weights on GPU

T = zeros(nf, nf); % Initialize T matrix
QI = zeros(nf, nf); % Initialize QI matrix

nv2 = gpuArray(single(repmat(nv, 1, 1, gq_number))); % Replicate normals
centroid = gpuArray(single(repmat(centroid, 1, 1, gq_number))); % Replicate centroids

for i = 1:nf
    % Define distance function
    r_square = @(x) sum((x - centroid(i,:,:)).^2, 2);

    % T matrix kernel function
    fT = @(x) 1./(4*pi) * r_square(x).^-0.5;

    % Calculate T matrix entries
    T(i,:) = gather(sum(W .* fT(Xq), [2, 3]));

    % QI matrix kernel function
    fQ = @(x) 1/(-4*pi) .* r_square(x).^(-3/2) .* sum((x - centroid(i,:,:)) .* nv2, 2);

    % Calculate QI matrix entries
    QI(i,:) = gather(sum(W .* fQ(Xq), [2, 3]));

    % Handle singularity for diagonal entries
    T(i,i) = singular_intergral(mesh.vertices(mesh.faces(i,1),:), ...
                                mesh.vertices(mesh.faces(i,2),:), ...
                                mesh.vertices(mesh.faces(i,3),:));
    QI(i,i) = 1/2;
end

% Solve the system to compute D2N
D2N = T \ QI;

fprintf('D2N calculation completed.\n');
end
