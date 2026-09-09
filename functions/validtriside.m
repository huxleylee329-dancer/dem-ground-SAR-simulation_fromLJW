function [TRI_idx_new,x_new,y_new,ref_PS_idx_new,valid_TRI_idx_new_log,valid_PS_idx_new_log] = validtriside(TRI_idx,x,y,ref_PS_idx,valid_TRI_idx)
%x为PS点的行索引，y为PS点的列索引，
%% 对PS点和PS网络进行重新排布,计算包含点第ref_PS_idx个PS点在内的最大连通三角网络（仅包含valid_PS_idx指示的PS点）
    if valid_TRI_idx(ref_PS_idx) == 0
        error('The reference PS must be valid!');
    end
    
    % 寻找原PS网络中包含第ref_PS_idx个PS点的最大连通网络，只有这个网络中的边才是有效边，且只有这些边的端点是有效的PS点
%     valid_TRI_idx_new_log = contriside(TRI_idx,ref_PS_idx) & valid_TRI_idx;
    valid_TRI_idx_new_log = contriside(TRI_idx,ref_PS_idx,valid_TRI_idx);
%     [~,~,ic] = unique(TRI_idx(valid_TRI_idx_new_log,:),'rows');
%     ind = (ic ~= 1);
    
    % 根据valid_TRI_idx_new_log，剔除无效的PS点
    valid_PS_idx_new = sort(unique(reshape(TRI_idx(valid_TRI_idx_new_log,:),[],1)),'ascend'); % 只有这些索引号的PS点是有效的
    
    valid_PS_idx_new_log = false(size(x));
    valid_PS_idx_new_log(valid_PS_idx_new) = 1;
    
    x_new = x(valid_PS_idx_new_log);
    y_new = y(valid_PS_idx_new_log);
    ref_PS_idx_new = find(valid_PS_idx_new == ref_PS_idx);

    ps_count_new = numel(x_new);
    
    TRI_idx_new_tmp = TRI_idx(valid_TRI_idx_new_log,:);
    TRI_idx_new = zeros(size(TRI_idx_new_tmp));
    
    for i = 1:ps_count_new
        TRI_idx_new(TRI_idx_new_tmp == valid_PS_idx_new(i)) = i;
    end
end