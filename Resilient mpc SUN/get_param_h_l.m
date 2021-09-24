function [h1,h2,lambda1,lambda2,ok] = get_param_h_l(A,B,K,P,N,M,coeff2)
%This script compute the parameters h and lambda as describe ine the report page 34
lambda1 = abs(max(eig(A)));
if coeff2
    lambda2 = 2*abs(max(eig(A+B*K)));
else
    lambda2 = abs(max(eig(A+B*K)));
end

norm2_Ak = []; lambda1k = [];
norm2_ABKk = []; lambda2k = [];

T = chol(P);
iT=inv(T);
for k=1:M
    norm2_Ak = [ norm2_Ak norm(T*A^k*iT,2)];
    lambda1k = [lambda1k lambda1^k];
end

for k=1:N-M
    norm2_ABKk = [ norm2_ABKk norm(T*(A+B*K)^k*iT,2)];
    lambda2k = [lambda2k lambda2^k];      
end

fact_h1 = norm2_Ak./lambda1k; fact_h2 = norm2_ABKk./lambda2k;

h1 = max(fact_h1);
h2 = max(fact_h2);


ok = f(M,N,h1,h2,lambda1,lambda2) < 1;

end