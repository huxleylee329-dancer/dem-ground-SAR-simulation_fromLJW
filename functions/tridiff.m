function PS_val_diff = tridiff(PS_val,TRI_idx)
    L = size(TRI_idx,1);
    PS_val_diff = zeros(L,1);
    
    for i = 1:size(TRI_idx,1)
        PS_val_diff(i) = PS_val(TRI_idx(i,1)) - PS_val(TRI_idx(i,2));
    end
end