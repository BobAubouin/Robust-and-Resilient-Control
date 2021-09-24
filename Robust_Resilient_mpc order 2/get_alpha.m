function [alpha,ok] = get_alpha(A,B,D,K,N,M,Zw,H1)
%This function find the smallest parameter alpha such that the ellipsoide
%define with matrix H1 and alpha is a muRPI. Moreover the ok flag is 1 when
%alpha exist and 0 else.
[h1,h2,lambda1,lambda2,ok] = get_param_h_l(A,B,K,H1,N,M,1);
if ~ok
    alpha = Inf;
else
    s = [1 1 -1 -1;
         -1 1 -1 1];
     W = [];
     vertices_zW = vertices(zonotope([Zw.c,Zw.R]));
     for i=1:size(vertices_zW,2)
         W(i) = sqrt((D*(vertices_zW(:,i)))'*H1*(D*(vertices_zW(:,i))));
     end
    wbar = max(W);
    alpha = (g(M,N,h1,h2,lambda1,lambda2)*wbar/(1-f(M,N,h1,h2,lambda1,lambda2)))^2;
end

end