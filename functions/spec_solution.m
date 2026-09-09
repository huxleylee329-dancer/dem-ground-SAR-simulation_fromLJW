function x = spec_solution(A,b)
    N = length(b);
    [U,S,V] = svd(A);
    right_part = U'*b;
    
    r = rank(A);
    S_diag = diag(S(1:r,1:r));
    S = diag(1./S_diag);
    
    x = V*[S*right_part(1:r);zeros(N-r,1)];
end