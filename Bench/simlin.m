function [Y,t,X] = simlin(A,Bu,Bv,C,T,u,v,x0)

kmax = length(T);
n = size(A,1);
X = zeros(n,kmax+1);
X(:,1) = x0;
for k = 1:kmax
    X(:,k+1) = A*X(:,k)+Bu*u(:,k)+Bv*v(:,k);
end

Y = C*X;
t = [2*T(1)-T(2),T];
end
