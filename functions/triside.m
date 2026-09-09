function L = triside(TRI)
N = size(TRI,1); % Number of triangles
idx = zeros(3*N,2);
TRI_sort = sort(TRI,2,'ascend');

for i = 1:N
    idx(3*i-2,:) = [TRI_sort(i,1),TRI_sort(i,2)]; % 1st side of the i-th triangle
    idx(3*i-1,:) = [TRI_sort(i,2),TRI_sort(i,3)]; % 2nd side of the i-th triangle
    idx(3*i,:)   = [TRI_sort(i,1),TRI_sort(i,3)]; % 3rd side of the i-th trianlge
end

L = unique(idx,'rows'); % Eliminate the same side and sort the rows  返回唯一的行数
end