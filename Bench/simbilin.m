function [Y,t,X] = simbilin(A,Bu,Bv,Bxu,Bvu,C,T,u,v,x0)

kmax = length(T);
n = size(A,1);
X = zeros(n,kmax+1);
X(:,1) = x0;
for k = 1:kmax
    int = 0; %variable intermédiaire
    for i=1:size(Bu,2)
        int = int + (Bxu(:,:,i)*X(:,k)+Bvu(:,:,i)*v(:,k))*u(i,k);
    end
    X(:,k+1) = A*X(:,k)+Bu*u(:,k)+Bv*v(:,k)+int;
end

Y = C*X;
t = [2*T(1)-T(2),T];
end
