function [test_ok,B] = TestMRPI(M,N,A,ABK,DW,RPI,H,B,Z0,nstep)
%This function test if a RPI is a mu-RPI ofr a number of attack M over the
%horizon prediction N. Moreovver it return the bound of the trajectory
%reach for each scenario of attack. We cosider the dynamic:
% x+ = (A + nu BK)x + Dw
% 
% where nu is the variable which models the attack (nu=0 when an attack occurs, nu=1 else)
% 
% Input:  - M: Number of attack
%         - N: Horizon prediction
%         - A: dynamic matrix of the system without feedback (when an attack oocurs)
%         - ABK: dynamic matrix (A+BK) of the system with feedback (without attack)
%         - DW: Set for the disturbances
%         - RPI : RPI set for the dynamic x+ = (A + BK)x + Dw
%         - Z0 : argument used for recursivity 
%         - nstep : argument used for recursivity 
% 
% Output: - test_ok: boolean = 1 if the condition is vberified, =0 else
%         - P_traj: extrem point of each set reachable for the scenarios
%         
% Bob Aubouin 30/09/2021


if nargin==8
    Z0 = RPI;% à l'instant 0 le système est dans le RPI
    nstep=1; % première attaque à l'instant k=1
end

if M==0 %Plus d'attaque à tester, on arrête la récursivité
    test_ok = 1;
else
    test_ok = 0;

    for k=1:nstep %Instant de l'attaque à tester
        step = Inf;
        Z = Z0; %On repart toujours du set atteint à l'attaque précédente
        for i=1:N
            if i==k %l'attaque arrive à l'instant k
                Z = A*Z + DW;
                Zint = Z;
                %Zint = reduce(Z,'combastel',10);
            else
                Z = (ABK)*Z + DW;
            end
            for j=1:length(B)-1
                Hj = H(:,j);
                B(j) = max(B(j), Hj'*Z.Z(:,1) + sum(abs(Hj'*Z.Z(:,2:end))));
            end
            B(end) = B(end)+1;
            if in(RPI,Z) % on vérifie si on est dans le RPI
                step = i; %Si oui on note l'étape et on arrête la boucle
                break
            end
        end
        if step<=N %Si on a réussi à atteindre le RPI en moins de N étape on peut tester les autre possibilité d'attaque
            %disp([M,N,step])
            [test_int,B_int] = TestMRPI(M-1,N-k,A,ABK,DW,RPI,H,B,Zint,step-1); %Partie Recursive, on test de nouveaux scénario avec 
            if ~test_int %Si on a pas atteint le RPI on s'arrête
                return
            end
            B = max(B,B_int);
        else %Si on a pas atteint le RPI on s'arrête
            return
        end      
    end
    test_ok = 1;
end