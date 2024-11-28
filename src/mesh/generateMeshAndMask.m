function [mesh, mask] = generateMeshAndMask(radbound, Mask, voxel_size, isiso)
% generateMeshAndMask
% Generate a mesh and a mask inside the mesh based on the iso2mesh toolbox
% (available at http://iso2mesh.sf.net). This function is only suitable for
% data with isotropic voxel size or anisotropic data where the slice resolution
% is approximately 1/3 of the in-plane resolution (e.g., [0.625, 0.625, 2]).
%
% INPUT:
%   radbound    - Radius boundary parameter for surface generation.
%   Mask        - 3D binary matrix representing the initial mask.
%   voxel_size  - Array specifying the voxel size in each direction [dx, dy, dz].
%   isiso       - Logical flag (1 or 0) indicating if the data has isotropic voxel size.
%
% OUTPUT:
%   mesh        - Struct containing:
%                   - mesh.vertices: Vertex coordinates (N x 3 matrix)
%                   - mesh.faces: Triangular mesh faces (M x 3 matrix)
%   mask        - 3D binary matrix indicating the mask inside the mesh
%
%   Created by Haodong zhong (zhonghaodonghy@outlook.com) on 2022.06
%   Last modified by Haodong zhong on 2024.11.25

matrix_size = size(Mask);
opt.radbound = radbound;%2.8
opt.distbound = 1; % Distance boundary parameter
inslice = voxel_size(1);
if isiso
    % Set parameters for surface generation
    Mask_model = Mask;
    Mask_model = SMV(Mask_model,size(Mask_model), [inslice,inslice,inslice], 1)>0.001;
    % Generate initial triangular surface from the mask
    [no, el, ~, ~] = vol2surf(double(Mask_model), 1:matrix_size(1), ...
        1:matrix_size(2), 1:matrix_size(3), opt, 0, 'cgalsurf', 1);
else
    Mask_model = zeros([matrix_size(1),matrix_size(2),matrix_size(3)*3]);
    Mask_model(:,:,1:3:end) =Mask;
    Mask_model(:,:,2:3:end) =Mask;
    Mask_model(:,:,3:3:end)=Mask;
    Mask_model = double(Mask_model);
    inslice = voxel_size(1);
    Mask_model = SMV(Mask_model,size(Mask_model), [inslice,inslice,inslice], 2)>0.001;
    Mask_model = SMV(Mask_model,size(Mask_model), [inslice,inslice,inslice], 1)>0.999;
    [no,el, ~, ~]=vol2surf(Mask_model,1:matrix_size(1),1:matrix_size(2),1:matrix_size(3)*3,opt,1,'cgalsurf',1);
    no(:,3)=no(:,3)./3+2/3;
    %Mask_model = Mask;
end
% Compute vertex coordinates
vertices =(no+[-(matrix_size(1)-1)/2-1,-(matrix_size(2)-1)/2-1,-(matrix_size(3)-1)/2-1]).*voxel_size';
faces = el(:,1:3);
% Compute the area of each triangle
area_tri = area(vertices(faces(:, 1), :), vertices(faces(:, 2), :), vertices(faces(:, 3), :));
% Remove degenerate triangles based on area threshold
faces(area_tri < 10^(-2), :) = [];

% Adjust the mesh to ensure normals are properly oriented
vertices = vertices+inslice/2;
mesh.vertices = double(vertices);
mesh.faces = faces;
mesh = adjust_mesh(mesh);
mesh = adjustMeshOrientation(mesh, Mask, voxel_size);
%% Check points inside the mesh
%maskb = Mask-(SMV(Mask,matrix_size, voxel_size, 4)>0.999);
%maskb = maskPointsInsideMesh(logical(maskb), voxel_size, mesh);
%in = inpolyhedron(mesh, X, Y, Z);
%in2 = double(in);
%mask0 = permute(in2, [2, 1, 3]);
%mask = maskb+(SMV(Mask_model,matrix_size, voxel_size, 4)>0.999);
%mask = maskb;
%mask = mask+(SMV(Mask_model,matrix_size, voxel_size, 4)>0.999);
mask = logical(Mask);

end
