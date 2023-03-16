function [alpha,Xf_e,Xf_z] = TerminalSet_P(K,P,Xu,Xx,N,Nz)
%This function calcul the terminal set Xf, the bigger ellipsoidal set under
%the form %Xf = {x | x'Px <alpha} such that Xf includes in Xx and K Xf includes in Xu

%Inputs :   - K: Feedback matrices
%           - P: matrice solution of the Liapunov equation P - (A-BK)'P(A-BK) = Q +
%           K'Rk describing the ellispoid
%           - Xu : constraint set for the command input
%           - Xx : constrain set for the states
%           - N: maximal numer of iteration for the dichotomi algorithm
%           (default N=10)
%           - Nz: Number of segment for the zonotopique approximation of
%           the final ellispoid (default Nz=20)
%Outputs:   - alpha: final alpha describing the ellipsoid
%           - Xf_e: final ellipsoid set
%           - Xf_z: zonotopique approximation of the final ellispoid set

if nargin == 4
    N = 10;
    Nz = 20;
elseif nargin == 5
    Nz = 20;
end

E0 = ellipsoid(inv(P));

%Initialisation
AlphaX = sqrt(normP(Xx)/norm(E0));
AlphaU = sqrt(normP(Xu)/norm(K*E0));

if AlphaX < AlphaU
    alphaplus = AlphaX;
else
    alphaplus = AlphaU;
end
alphamoins = 0;
E_alpha = ellipsoid(inv(P)*(alphaplus-alphamoins)/2);

i = 0; %compteur
alpha = (alphaplus-alphamoins)/2;
while ~in(Xx,E_alpha) && i<N   %boucle par dichotomie
    E_alpha = ellipsoid(inv(P)*alpha);  
    if in(Xx,E_alpha) && in(Xu,K*E_alpha) 
        alphamoins = alpha;
    else
        alphaplus = alpha;
    end
    alpha = (alphaplus-alphamoins)/2;
    i = i+1;
end

Xf_e = ellipsoid(inv(P)*alpha); 
Xf_z = zonotope(Xf_e,Nz,'i:norm'); %Zonotopique approximation

end



    
    
