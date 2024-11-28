function restored_matrix = u_map(u,mask)
% 原始三维矩阵的大小
original_size = size(mask);

% 创建一个与原始矩阵大小相同的全零矩阵
restored_matrix = zeros(original_size);

% 根据 mask 非零位置将列向量 f 填回矩阵
restored_matrix(logical(mask)) = u;

end