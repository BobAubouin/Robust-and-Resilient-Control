function [DoS] = create_attack(N,M)
% This function returns a stack of vector representing all the possiblity
% for M attack during an interval of size N. (0 = attack, 1= normal)
DoS = zeros(N,nchoosek(N,M));
id = nchoosek(1:N,M);

for i=1:size(id,1)
   DoS(id(i,:),i) = 1;    
end
DoS = 1 - DoS;
end

