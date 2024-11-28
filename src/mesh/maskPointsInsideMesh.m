function in_mask = maskPointsInsideMesh(mask, voxel_size, mesh)
    %MASKPOINTSINSIDEMESH Determines which voxels in a mask are inside a mesh.
    %
    %   IN_MASK = MASKPOINTSINSIDEMESH(MASK, VOXEL_SIZE, MESH) returns a logical
    %   array IN_MASK of the same size as MASK, where IN_MASK(i,j,k) is true if the
    %   voxel at position (i,j,k) is inside the mesh defined by MESH.
    %
    %   Inputs:
    %       - mask: a 3D logical array representing the region of interest.
    %       - voxel_size: a 1x3 vector [vx, vy, vz] specifying the size of each voxel.
    %       - mesh: a structure with fields 'vertices' and 'faces' defining the mesh.
    %
    %   Output:
    %       - in_mask: a 3D logical array of the same size as mask, where true values
    %         indicate that the voxel is inside the mesh.

    % 确保 mask 为逻辑类型
    mask = logical(mask);

    % 获取 mask 的尺寸
    [Nx, Ny, Nz] = size(mask); % 注意，这里是 [X，Y，Z]

    % 生成网格坐标
    gridx0 = -(Nx-1)/2 : (Nx-1)/2;
    gridy0 = -(Ny-1)/2 : (Ny-1)/2;
    gridz0 = -(Nz-1)/2 : (Nz-1)/2;

    % 缩放网格坐标以匹配物理尺寸
    gridx = gridx0 * voxel_size(1);
    gridy = gridy0 * voxel_size(2);
    gridz = gridz0 * voxel_size(3);

    % 使用 ndgrid 生成坐标网格，注意坐标轴的顺序
    [X, Y, Z] = ndgrid(gridx, gridy, gridz); % 调整为 [X, Y, Z]

    % 提取 mask 中为 true 的坐标点
    X_masked = X(mask);
    Y_masked = Y(mask);
    Z_masked = Z(mask);

    % 组合成查询点矩阵
    query_points = [X_masked(:), Y_masked(:), Z_masked(:)];

    % 判断点是否在 mesh 内部
    in = inpolyhedron(mesh, query_points);

    % 初始化输出的 in_mask
    in_mask = false(size(mask));

    % 将结果映射回对应的体素位置
    in_mask(mask) = in;

%     % 绘制 mesh
%     % Plot the mesh
%     patch('Faces', mesh.faces, 'Vertices', mesh.vertices, ...
%           'FaceAlpha', 0.2, 'LineStyle', ':', 'EdgeColor', 'blue');
%     daspect([1 1 1]);
%     view(3);
%     axis tight;
%     set(gca, 'YTick', -1:1:1);
%     camlight;
%     hold on;
% 
% % Extract the coordinates of points where mask is 1
% pointsToPlot = [X(mask), Y(mask), Z(mask)]; % Use mask directly to filter X, Y, Z
% 
% % Plot the points
% scatter3(pointsToPlot(:,1), pointsToPlot(:,2), pointsToPlot(:,3), 'g.');




end

