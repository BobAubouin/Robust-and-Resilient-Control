%% version ligne de commandes
A = [1.1 2; 0 0.95];
B = [0; 0.0787];
C = [-1 1];
D = 0;
Ts = 1;
sys = ss(A,B,C,D,Ts);
x0 = [0.5;-0.5]; % initial states at [0.5 -0.5]

Qy = 1;
R1 = 1;
R2 =0.5;
Ref = 1;
K_lqr = lqry(sys,Qy,R1);

M = [C*A;C*A^2;C*A^3;C*A^4];
CONV = [C*B 0  0 0;...
        C*A*B C*B 0 0;...
        C*A^2*B C*A*B C*B 0;...
        C*A^3*B C*A^2*B C*A*B C*B];

Q = Qy;
Q_bar = dlyap((A-B*K_lqr)', Q+K_lqr'*R1*K_lqr);
Q_hat = blkdiag(Q,Q,Q,Q);
R_hat = blkdiag(R1,R1,R1,R1);
REF = [Ref;Ref;Ref;Ref];

U2DU = [1 0 0 0;
        -1 1 0 0;
        0 -1 1 0;
        0  0 -1 1;];
 U2U0 = [1 0 0 0];
    
DR = U2DU'*blkdiag(R2,R2,R2,R2)*U2DU;

H = CONV'*Q_hat*CONV + R_hat+DR;



L = chol(H,'lower');
Linv = L\eye(size(H,1));


Ac = [1 0 0 0;...
      -1 0 0 0;...
       0 1 0 0;...
       0 -1 0 0;...
       0 0 1 0;...
       0 0 -1 0;...
       0 0 0 1;...
       0 0 0 -1];
b0 = [1;1;1;1;1;1;1;1];


x = x0;
iA = false(size(b0));

% create options for the solver, and specify non-hessian first input
opt = mpcActiveSetOptions;
opt.IntegrityChecks = false;
opt.UseHessianAsInput = false;
x = x0;
t_constrained = 0:40;
u = [0 0 0 0];
for ct = t_constrained
    yMPC(ct+1) = C*x;
    F = CONV'*Q_hat*M*x - (REF'*Q_hat*CONV)' - (u(1)'*R2*U2U0)';
    [u,status,iA] = mpcActiveSetSolver(Linv,F,Ac,b0,[],zeros(0,1),iA,opt);
    uMPC(ct+1) = u(1);
    x = A*x+B*uMPC(ct+1);
end



%% mpc matlab

iMV = 1;
iMO = 1;

Plant = setmpcsignals(sys,'MV',iMV,'MO',iMO);


MPCobj = mpc(Plant,Ts);

MPCobj.PredictionHorizon = 4; % 4h
MPCobj.ControlHorizon = 4;

MPCobj.Weights.ManipulatedVariables = sqrt(R1);
MPCobj.Weights.OutputVariables = sqrt(Qy);
MPCobj.Weights.ManipulatedVariablesRate = sqrt(R2);


MPCobj.ManipulatedVariables(1).Min = -1;
MPCobj.ManipulatedVariables(1).Max = 1;
xmpc =  mpcstate(MPCobj);
x = x0;
t_constrained = 0:40;

for ct = t_constrained
    yMPC2(ct+1) = C*x;
    xmpc.Plant = x;
    u = mpcmove(MPCobj,xmpc,yMPC2(ct+1),Ref);
    uMPC2(ct+1) = u;
    x = A*x+B*uMPC2(ct+1);  
end

figure
subplot(2,1,1)
plot(t_constrained,uMPC)
hold on
plot(t_constrained,uMPC2)
xlabel('time')
ylabel('u')
subplot(2,1,2)
plot(t_constrained,yMPC)
hold on
plot(t_constrained,yMPC2)
xlabel('time')
ylabel('y')
legend('diy mpc','matlab mpc')

