function Xn = InverseReach(Xf,A,B,Xu,Xx,N,nbDoS)
%%Version with the outter intersection with the state constraint
%
%
%Calcul the Initial feasible set, given a system dynamic and a terminal
%set. Also include the possibility to add DoS at the beginning of the
%command sequence
%Inputs:    - Xf: Terminal set (define with CORA)
%           - A,B: system matrices such that: x_{k+1} = A x_k + Bu
%           - Xu: constrain set for the control input
%           - Xx: constrain set for the state
%           - N: prediction horizon of the MPC
%           - Number of DoS attacks
%Outputs:   - Xn : Terminal feasible set (as a CORA set).
%
% Bob Aubouin, 23/06/2021

Xn = Xf;
for i=1:N
    if i<=N-nbDoS
        Xn = inv(A)*(Xn+(-1*B*Xu));
    else
      Xn = inv(A)*Xn;  
    end
    Xn = Xn & Xx;
end

end

