function PS_val = lsqadjust(PS_Diff_val,TRI_idx,weight,Abs_ref_idx,Abs_ref_val)
%最小二乘平差法
    L = size(TRI_idx,1); % 三角网格的边数
    K = max(TRI_idx(:));
    
    B_PS = zeros(L,K);
    B_PS(sub2ind([L,K],(1:L)',TRI_idx(:,1))) = 1;
    B_PS(sub2ind([L,K],(1:L)',TRI_idx(:,2))) = -1;
    
    Coef_A = B_PS'*diag(weight.^2)*B_PS;
    Coef_b = B_PS'*diag(weight.^2)*PS_Diff_val;
    
    sps = spec_solution(Coef_A,Coef_b);
    PS_val = sps + (Abs_ref_val - sps(Abs_ref_idx))*ones(K,1);
end