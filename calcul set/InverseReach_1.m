function Xn = InverseReach_1(Xf,A,B,Xu,N,nbDoS)
%Version without the intersection with the state constraint
%
%Calcul the Initial feasible set, given a system dynamic and a terminal
%set. Also include the possibility to add DoS at the beginning of the
%command sequence
%Inputs:    - Xf: Terminal set (define with CORA)
%           - A,B: system matrices such that: x_{k+1} = A x_k + Bu
%           - Xu: constrain set for the control input
%           - N: prediction horizon of the MPC
%           - Number of DoS attacks
%Outputs:   - Xn : Terminal feasible set (as a CORA set).
V = vertices(Xf);
Xn = mptPolytope(Polyhedron('V',V'));
iA = inv(A);
for i=1:N
    if i<=N-nbDoS
        disp(Xn)
        Xn = iA*(Xn+(-1*B*Xu));
    else
      Xn = iA*Xn;  
    end
end

end

