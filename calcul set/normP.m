function n = normP(Poly)
%Compute the euclidean norm of a Polytope
V = vertices(Poly);
n = 0;
for i=1:size(V,2)
    n = max(n,norm(V(:,i)));
end
end