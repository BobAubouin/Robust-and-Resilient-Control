function [u,new_data,Xk] = MPC_maison(mpc,dY,Xk_1,Ref,Vm,u_prec,M,prev_Data,Zono)
%% system restriction data
Q = diag(mpc.Weights.OutputVariables.^2);
R = diag(mpc.Weights.ManipulatedVariables.^2);
dR = diag(mpc.Weights.ManipulatedVariablesRate.^2);


A = mpc.Model.Plant.A;
B = mpc.Model.Plant.B;
C = mpc.Model.Plant.C;
Bv = B(:,mpc.Model.Plant.InputGroup.Measured);
Bu = B(:,mpc.Model.Plant.InputGroup.Manipulated);
Bd = B(:,mpc.Model.Plant.InputGroup.Unmeasured);

N = mpc.PredictionHorizon;
n = size(A,1);
nu = size(Bu,2);
nv = size(Bv,2);
nc = size(C,1);

Vm = reshape(Vm,[size(Vm,1),size(Vm,3)]);
Vm = Vm(:,1:min(N+1,size(Vm,2)));
if size(Vm,2)==N+1
   Vm = reshape(Vm,[(N+1)*nv,1]);
else
   Vm = [Vm,repmat(Vm(:,end),1,N+1-size(Vm,2))];
   Vm = reshape(Vm,[(N+1)*nv,1]);
end

%% Observer

%Xk = A*Xk_1 + Bu*u_prec + Bv*Vm(1:nv) +M*(dY-C*Xk_1);%
Xk = Xk_1;

Vm = Vm(nv+1:end);
%% define cost
if ~isfield(prev_Data,'Linv')
    R_hat = [];
    Q_hat = [];
    dR_hat = [];
    dRmove = eye(N*nu)+diag(-1*ones((N-1)*nu,1),-nu);

    for i=1:N
        R_hat = blkdiag(R_hat,R);
        dR_hat = blkdiag(dR_hat,dR);
        Q_hat = blkdiag(Q_hat,Q);
    end
    dR_hat = dRmove'*dR_hat*dRmove;

    Cu = [];
    Cv = [];
    for i=1:N
        Cu = blkdiag(Cu,C*Bu);
        Cv = blkdiag(Cv,C*Bv);
    end
    M = zeros(nc*N,n);
    M(1:nc,:) = C*A;

    for i=1:N-1
        tempu = [];
        tempv = [];
        for k=1:N-i
            tempu = blkdiag(tempu,C*A^i*Bu);
            tempv = blkdiag(tempv,C*A^i*Bv);
        end
        Cu(i*nc+1:end,1:(N-i)*nu) = Cu(i*nc+1:end,1:(N-i)*nu) + tempu;
        Cv(i*nc+1:end,1:(N-i)*nv) = Cv(i*nc+1:end,1:(N-i)*nv) + tempv;
        M(i*nc+1:(i+1)*nc,:) = C*A^(i+1);
    end
    
    H = Cu'*Q_hat*Cu +  R_hat + dR_hat;
    L = chol(H,'lower');
    Linv = L\eye(size(H,1));
    
    iA = false(2*N*nc+2*N*nu,1);
    
    prev_Data.Linv = Linv;
    prev_Data.H = H;
    prev_Data.M = M;
    prev_Data.Cv = Cv;
    prev_Data.Cu = Cu;
    prev_Data.Q_hat = Q_hat;
    prev_Data.iA = iA;
    
else
    Linv = prev_Data.Linv;
    H = prev_Data.H;
    M = prev_Data.M;
    Cv = prev_Data.Cv;
    Cu = prev_Data.Cu;
    Q_hat = prev_Data.Q_hat;
    iA = prev_Data.iA;
end

Rf = zeros(nc*N,1);
for i=0:N-1
    Rf(i*nc+1:(i+1)*nc) = Ref;
end

Beta = M*Xk+Cv*Vm;

U2u0 = zeros(nu,N*nu);
U2u0(:,1:nu) = eye(nu);


F = Cu'*Q_hat*Beta - Cu'*Q_hat*Rf - U2u0'*dR*u_prec;

%% Zono Constraints

Fx = zeros(2*N*nc+2*N*nu,N*nc);
Fu = zeros(2*N*nc+2*N*nu,N*nu);
Bineq = zeros(2*N*nc+2*N*nu,1);

dFx = zeros(2*nc,nc);
dFu = zeros(2*nu,nu);
dBx = zeros(nc,1);
dBu = zeros(nu,1);


for i=1:nc
    dFx(i*2-1:2*i,i) = [1;-1];
    dBx(i*2-1) = mpc.OutputVariables(i).Max;
    dBx(i*2) = -mpc.OutputVariables(i).Min;
end

for i=1:nu
    dFu(i*2-1:2*i,i) = [1;-1];
    dBu(i*2-1) = mpc.ManipulatedVariables(i).Max;
    dBu(i*2) = -mpc.ManipulatedVariables(i).Min;
end

for i=0:N-1
    Fx(2*i*nc+1:2*(i+1)*nc,i*nc+1:(i+1)*nc) = dFx;
    Fu(2*N*nc+2*i*nu+1:2*N*nc+2*(i+1)*nu,i*nu+1:(i+1)*nu) = dFu;
    Bineq(2*i*nc+1:2*(i+1)*nc) = dBx;
    Bineq(2*N*nc+2*i*nu+1:2*N*nc+2*(i+1)*nu) = dBu;
end

Ac = Fx*Cu+Fu;
Bc = Bineq-Fx*Beta;

% Ac = Ac(2*N*nc+1:end,:);
% Bc = Bc(2*N*nc+1:end,:);
opt = mpcActiveSetOptions;
opt.UseHessianAsInput = false;

%[U,~,iA] = quadprog(Linv,F,Ac,Bc);
[U,~,iA] = mpcActiveSetSolver(Linv,F,Ac,Bc,zeros(0,N*nu),zeros(0,1),iA,opt);
prev_Data.iA = iA;

u = U(1:nu);
q1=0.2; %discrétisation de la commande des volets
u(1:2)= q1 * round(u(1:2)/q1);
q2=5; %discrétisation de la commande des radiateurs
%u(3:5)= q2 * round(u(3:5)/q2);

new_data = prev_Data;
end