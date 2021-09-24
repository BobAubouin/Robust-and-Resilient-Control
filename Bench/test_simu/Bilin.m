function [X_k1,Y] = Bilin(X_k,u,v,A_k,Bu_k,Bv_k,Bxu_k,Bvu_k,C_k)
%Simulation d'un pas de temps d'un modèle bilinear
int = zeros(size(X_k,1),1); %variable intermédiaire
for i=1:size(Bu_k,2)
    int = int + (Bxu_k(:,:,i)*X_k + Bvu_k(:,:,i)*v)*u(i);
end
X_k1 = A_k*X_k + Bu_k*u + Bv_k*v + int;
Y = C_k*X_k1;
end
