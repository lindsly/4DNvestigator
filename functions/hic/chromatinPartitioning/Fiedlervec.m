function [FN,Fvec] = Fiedlervec(A)
    A = A - diag(diag(A)); A = 0.5*(A + A.');
    [dummy, L_norm] = normalization_A(A); eps = 1e0;
    if length(A) > 1
        [Eig2cols Deig] = eigs(sparse(L_norm+eps*eye(size(A,1))), 2, 'SM');
        idx_2nd = find(diag(Deig)> ( eps + 1e-4) );
        if ~isempty(idx_2nd)
           Fvec = Eig2cols(:,idx_2nd);
           FN = Deig(idx_2nd,idx_2nd) - eps;
        else
            Fvec = []; 
            FN = Deig(end,end) - eps;
        end
    else
        FN = 0;
        Fvec = [];
    end
    
    if abs(FN) < 1e-10
        FN = 0;
    end
    
end