function [A_norm, L_norm] = normalization_A(A)
         A = 0.5*(A+A.');
         A_zd = A - diag(diag(A));   %%% remove diagonal for conventional adj. matrix
         deg_tmp = sum(A_zd,2); %%% degree matrix
         deg_tmp_sqrtinv = 1./sqrt(deg_tmp); 
         diag_ones = ones(length(deg_tmp_sqrtinv),1);
         diag_ones(isinf(deg_tmp_sqrtinv))=0; diag_ones(isnan(deg_tmp_sqrtinv))=0;
         deg_tmp_sqrtinv(isinf(deg_tmp_sqrtinv))=0; deg_tmp_sqrtinv(isnan(deg_tmp_sqrtinv))=0;
%          A_norm =  diag(deg_tmp_sqrtinv) * A_zd * diag(deg_tmp_sqrtinv);
         A_norm =  (deg_tmp_sqrtinv*deg_tmp_sqrtinv.').* A_zd;
         L_norm = diag(diag_ones) - A_norm;
end