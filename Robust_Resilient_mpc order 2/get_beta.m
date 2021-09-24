function Beta = get_beta(A,B,D,K,N,M,Zw,H2,alpha_opt,H1)

[h1,h2,lambda1,lambda2] = get_param_h_l(A,B,K,H2,N,M,1);
s = [1 1 -1 -1;
     -1 1 -1 1];
 W = [];
vertices_zW = vertices(zonotope([Zw.c,Zw.R]));
for i=1:size(vertices_zW,2)
    W(i) = sqrt((D*(vertices_zW(:,i)))'*H1*(D*(vertices_zW(:,i))));
end
wbar = max(W);
L = [];
T1 = chol(H1);
T2 = chol(H2);
int = norm(T2/T1,2);
for i=1:N
    L = [L (f(M,i,h1,h2,lambda1,lambda2)*int*sqrt(alpha_opt)+g(M,i,h1,h2,lambda1,lambda2)*wbar)^2];
end
Beta = max(L);

end