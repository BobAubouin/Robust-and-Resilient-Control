function [u,new_data,Xk] = MPC_zono(mpc,Xk_1,Vm,u_prec,prev_Data,Zono)
%% system restriction data
Q = mpc.Q;
R = mpc.R;
dR = mpc.dR;


A = mpc.Plant.A;
B = mpc.Plant.B;
C = mpc.Plant.C;
Bd = B(:,mpc.iUD);
if ~isempty(mpc.iMD)
    if model==1
        Bv = [B(:,mpc.iMD)];%, Bd];
    else
        Bv = [Bd, B(:,mpc.iMD)];
    end    
else
    Bv = [];
end
Bu = B(:,mpc.iMV);



    
N = mpc.PredictionHorizon;
n = size(A,1);
nu = size(Bu,2);
nv = size(Bv,2);
nc = size(C,1);

if isempty(Bv)
    Bv = zeros(n,1);
end

if ~isempty(Vm)
    Vm = reshape(Vm,[size(Vm,1),size(Vm,3)]);
    Vm = Vm(:,1:min(N+1,size(Vm,2)));
    if size(Vm,2)==N+1
       Vm = reshape(Vm,[(N+1)*nv,1]);
    else
       Vm = [Vm,repmat(Vm(:,end),1,N+1-size(Vm,2))];
       Vm = reshape(Vm,[(N+1)*nv,1]);
    end
    Vm_pred = Vm(nv+1:end);
else
    Vm_pred = zeros(N*1,1);
    Vm = 0;
    nv = 1;
end

%% Observer

Xk = Xk_1;

%% Zono Constraints
c_u = repmat(Zono.cu,N,1);
c_x = repmat(Zono.cx, N-1 , 1);
c_x = [c_x ; Zono.cterm];
R_u = [];
R_x = [];


for i=1:N-1
    R_u = blkdiag(R_u,Zono.Ru);
    R_x = blkdiag(R_x,Zono.Rx);
end

R_u = blkdiag(R_u,Zono.Ru);
R_x = blkdiag(R_x,Zono.Rterm);

nsu = size(Zono.Ru,2)*N;
nsx = size(R_x,2);
S2Sx = [zeros(nsx,nsu),eye(nsx)];
S2Su = [eye(nsu),zeros(nsu,nsx)];
    
    
if ~isfield(prev_Data,'Linv')
    %% define cost
    R_hat = [];
    Q_hat = [];
    dR_hat = [];
    dRmove = eye(N*nu)+diag(-1*ones((N-1)*nu,1),-nu);

    for i=1:N-1
        R_hat = blkdiag(R_hat,R);
        dR_hat = blkdiag(dR_hat,dR);
        Q_hat = blkdiag(Q_hat,Q);
    end
    R_hat = blkdiag(R_hat,R);
    dR_hat = blkdiag(dR_hat,dR);
    Q_hat = blkdiag(Q_hat,mpc.P);
    
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
    
    H = S2Sx'*R_x'*Q_hat*R_x*S2Sx + S2Su'*R_u'*(R_hat + dR_hat)*R_u*S2Su;
    H = H + 0.00001*eye(size(H,1));
    
    %% zono constraints
    Aineq = zeros((nsu+nsx)*2,(nsu+nsx));
    Bineq = ones((nsu+nsx)*2,1);

    for i=1:(nsu+nsx)
        Aineq(2*i-1,i)= 1;
        Aineq(2*i,i)= -1;
    end
    
    Aeq = Cu*R_u * S2Su - R_x*S2Sx;
    
    

    
    L = chol(H,'lower');
    Linv = L\eye(size(H,1));
    iA = false((nsu+nsx)*2,1);
    
    prev_Data.Linv = Linv;
    prev_Data.H = H;
    prev_Data.Aineq = Aineq;
    prev_Data.Bineq = Bineq;
    prev_Data.Aeq = Aeq;
    prev_Data.iA = iA;
    prev_Data.M = M;
    prev_Data.Cv = Cv;
    prev_Data.Cu = Cu;
    prev_Data.Q_hat = Q_hat;
    prev_Data.R_hat = R_hat;
    prev_Data.dR_hat = dR_hat;
else
    Linv = prev_Data.Linv;
    H = prev_Data.H;
    Aineq = prev_Data.Aineq;
    Bineq = prev_Data.Bineq;
    Aeq = prev_Data.Aeq;
    iA = prev_Data.iA;
    M = prev_Data.M;
    Cv = prev_Data.Cv;
    Cu = prev_Data.Cu;
    Q_hat = prev_Data.Q_hat;
    R_hat = prev_Data.R_hat;
    dR_hat = prev_Data.dR_hat;
end

Beta = M*Xk+Cv*Vm_pred;
Beq = c_x - Beta - Cu*c_u;


U2u0 = zeros(nu,nu*N);
U2u0(:,1:nu) = eye(nu);

F = S2Sx'*R_x'*Q_hat*c_x + S2Su'*R_u'*(R_hat+dR_hat)*c_u - S2Su'*R_u'*U2u0'*dR*u_prec;

% Ac = Ac(2*N*nc+1:end,:);
% Bc = Bc(2*N*nc+1:end,:);
opt = mpcActiveSetOptions;
opt.UseHessianAsInput = false;

%[U,~,iA] = quadprog(Linv,F,Ac,Bc);
[S,info,iA] = mpcActiveSetSolver(Linv,F,Aineq,Bineq,Aeq,Beq,iA,opt);
prev_Data.iA = iA;
if info==-1 || info==-2
    disp("problem")
end

u = Zono.Ru*S(1:nsu/N) + Zono.cu;


new_data = prev_Data;

% Xplot = zeros(n,N+1);
% Xplot (:,1) = Xk;
% Sx = S2Sx*S;
% for i=2:N
%     %Xplot (:,i) = Zono.Rx_mpc*Sx(12*i-23:12*i-12) + Zono.cx_mpc;
%     Xplot (:,i) = Zono.Rx*Sx(2*i-3:2*i-2) + Zono.cx;
% end  
% Xplot (:,N+1) = Zono.Rterm*Sx(19:23) + Zono.cterm;
% plot(Xplot(1,:),Xplot(2,:),'o-')
% hold on
% plot(Xplot(1,1),Xplot(2,1),'ro')
% grid on
end