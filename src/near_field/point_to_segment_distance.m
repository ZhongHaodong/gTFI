function distances = point_to_segment_distance(P, A, B)
% 计算点 P 到线段 AB 的距离
% P 是 N x 3
% A 和 B 是 1 x 3 的向量

% 从 A 到 B 的向量
AB = B - A;

% 从 A 到 P 的向量
AP = bsxfun(@minus, P, A);

% AP 在 AB 上的投影
t = sum(bsxfun(@times, AP, AB), 2) / sum(AB.^2);

% 将 t 限制在 [0,1] 范围内
t = max(min(t, 1), 0);

% 线段上的最近点
closest = bsxfun(@plus, A, t * AB);

% 距离
distances = sqrt(sum((P - closest).^2, 2));
end