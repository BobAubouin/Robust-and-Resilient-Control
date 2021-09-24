function [u,new_data,Xn] = MPC_ROBUSTE_improved(mpc,Xk,Ref,Vm,u_prec,prev_Data)


%% system restriction data
Q = mpc.Q;
R = mpc.R;
dR = mpc.dR;

K = mpc.K;
P = mpc.P;

A = mpc.Plant.A;
B = mpc.Plant.B;
C = mpc.Plant.C;
if ~isempty(mpc.iMD)
    Bv = B(:,mpc.iMD);
else
    Bv = [];
end
Bu = B(:,mpc.iMV);
Bd = B(:,mpc.iUD);


    
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
    nv = 1;
end
%% Calcul de la référence:

if ~isfield(prev_Data,'u0')
    if ~any(Ref)
        x0 = zeros(n,1);
        u0 = zeros(nu,1);
    else
        cvx_begin
            variables x0(n,1) u0(nu,1) r(nc,1)
            minimize(u0'*R*u0 + (r-Ref)'*(Q*10)*(r-Ref))
            subject to
                [A-eye(n) Bu; C zeros(nc,nu)] * [x0;u0] == [zeros(n,1);r];
        cvx_end
    end
    prev_Data.u0 = u0;
    prev_Data.x0 = x0;
else
    u0 = prev_Data.u0;
    x0 = prev_Data.x0;
end




%% Zono Constraints
Zono = mpc.zono;
c_u = repmat(Zono.cu_mpc - u0 , N , 1);
c_x = repmat(Zono.cx_mpc - x0 , N-1 , 1);
c_x = [c_x ; Zono.cterm];
R_u = [];
R_x = [];

for i=1:N-1
    R_u = blkdiag(R_u,Zono.Ru_mpc);
    R_x = blkdiag(R_x,Zono.Rx_mpc);
end

R_u = blkdiag(R_u,Zono.Ru_mpc);
R_x = blkdiag(R_x,Zono.Rterm);

nsu = size(Zono.Ru_mpc,2)*N;
nsxn = size(Zono.Rterm,2);
nsx = size(Zono.Rx_mpc,2)*(N-1) + nsxn;
nsx0 = size(Zono.Rz,2);
nsu0 = size(Zono.Ru_mpc,2);


S2sx0 = [zeros(nsx0,nsu),eye(nsx0),zeros(nsx0,nsx)];
S2sx1 = [zeros(nsx,nsu+nsx0),eye(nsx)];
S2su = [eye(nsu),zeros(nsu,nsx+nsx0)];
S2su0 = [eye(nsu0),zeros(nsu0,nsu-nsu0+nsx0+nsx)];
    
if ~isfield(prev_Data,'Linv')    
    
    %% define cost
    R_hat = [];
    Q_hat = [];
    dR_hat = [];
    dRmove = eye(N*nu)+diag(-1*ones((N-1)*nu,1),-nu);
    dRmove = [dRmove;zeros(nu,(N-1)*nu),eye(nu)];
    
    for i=1:N-1
        R_hat = blkdiag(R_hat,R);
        dR_hat = blkdiag(dR_hat,dR);
        Q_hat = blkdiag(Q_hat,C'*Q*C);
    end
    R_hat = blkdiag(R_hat,R);
    dR_hat = blkdiag(dR_hat,dR);
    dR_hat = blkdiag(dR_hat,dR);% traille nu(N+1) pour lui !
    Q_hat = blkdiag(Q_hat,P); %Gain final pour la stabilité
    dR_hat = dRmove'*dR_hat*dRmove;

    Cu = [];
    Cv = [];
    for i=1:N
        Cu = blkdiag(Cu,Bu);
        Cv = blkdiag(Cv,Bv);
    end
    M = zeros(n*N,n);
    M(1:n,:) = A;

    for i=1:N-1
        tempu = [];       
        tempv = [];
        for k=1:N-i
            tempu = blkdiag(tempu,A^i*Bu);
            tempv = blkdiag(tempv,A^i*Bv);
        end
        Cu(i*n+1:end,1:(N-i)*nu) = Cu(i*n+1:end,1:(N-i)*nu) + tempu;
        Cv(i*n+1:end,1:(N-i)*nv) = Cv(i*n+1:end,1:(N-i)*nv) + tempv;
        M(i*n+1:(i+1)*n,:) = A^(i+1);
    end
    
    H = S2sx1'*R_x'*Q_hat*R_x*S2sx1 + S2sx0'*Zono.Rz'*C'*Q*C*Zono.Rz*S2sx0 + S2su'*R_u'*(R_hat+dR_hat)*R_u*S2su;
    H = 0.00000000001*eye(size(H,1)) + H;
    
    %% zono constraints
    Aineq = zeros((nsu+nsx+nsx0)*2,(nsu+nsx+nsx0));
    Bineq = ones((nsu+nsx+nsx0)*2,1);

    for i=1:(nsu+nsx+nsx0)
        Aineq(2*i-1,i)= 1;
        Aineq(2*i,i)= -1;
    end
    
    Aeq = R_x*S2sx1 + M*Zono.Rz*S2sx0 - Cu*R_u*S2su;
    
    

    
    L = chol(H,'lower');
    Linv = L\eye(size(H,1));
    iA = false((nsu+nsx+nsx0)*2,1);
    
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
    
    %% modele nominal 
    
    
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

Beq = Cv*Vm_pred - c_x + Cu*c_u + M*(Xk - Zono.cz - x0);

F = S2sx1'*R_x'*Q_hat*c_x + S2su'*R_u'*(R_hat+dR_hat)*c_u - S2sx0'*Zono.Rz'*C'*Q*C*(Xk - Zono.cz - x0);% - S2su0'*Zono.Ru_mpc'*dR*u_prec;

% Ac = Ac(2*N*nc+1:end,:);
% Bc = Bc(2*N*nc+1:end,:);  
opt = mpcActiveSetOptions;
opt.UseHessianAsInput = false;

%[U,~,iA] = quadprog(Linv,F,Ac,Bc);
[S,info,iA] = mpcActiveSetSolver(Linv,F,Aineq,Bineq,Aeq,Beq,iA,opt);
prev_Data.iA = iA;

if info==-1
    disp("problem")
end
umpc = Zono.Ru_mpc*S(1:nsu0) + Zono.cu_mpc;

prev_Data.u_prec = umpc;

%% robuste
u = prev_Data.u_prec - K*(Zono.Rz*S2sx0*S + Zono.cz);

Xn = Xk - Zono.Rz*S2sx0*S - Zono.cz ;

% q1=0.2; %discrétisation de la commande des volets
% u(1:2)= q1 * round(u(1:2)/q1);
% q2=5; %discrétisation de la commande des radiateurs
% u(3:5)= q2 * round(u(3:5)/q2);

%% plot
% Xplot = zeros(n,N+1);
% Xplot (:,1) = Xk - Zono.Rz*S2sx0*S - Zono.cz;
% Sx = S2sx1*S;
% for i=2:N+1
%     Xplot (:,i) = Zono.Rx_mpc*Sx(2*i-3:2*i-2) + Zono.cx_mpc;
% end  
% plot(Xplot(1,:),Xplot(2,:))
% hold on
% plot(Xplot(1,1),Xplot(2,1),'ro')
% grid on



new_data = prev_Data;
end