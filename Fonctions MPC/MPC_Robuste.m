function [U,new_data,Xk] = MPC_Robuste(mpc,Xk_1,Ref,Vm,u_prec,prev_Data,model)

%% A faire :

% Calculer le coup d'arrivée dans notre cas : Fait à moitié
% Ajouter la simulation du modèle nominal : OK
% Faire la somme de la commande proportionnel et mpc: OK
% Faire le calcul du set invariant pour Xf : Pas fait
% Contraindre Xn dans le set Xf : Pas fait

%% system restriction data
Q = mpc.Q;
R = mpc.R;
dR = mpc.dR;

K = mpc.K;
P = mpc.P;


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
    prev_Data.Vm = 0;
    Vm = 0;
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


%% Observer

Xk = Xk_1;
    

%% Zono Constraints
Zono = mpc.zono;
c_u = repmat(Zono.cu_mpc-u0,N,1);
c_x = repmat(Zono.cx_mpc-x0, N-1 , 1);
c_x = [c_x ; Zono.cterm];
R_u = [];
R_x = [];

for i=1:N-1
    R_u = blkdiag(R_u,Zono.Ru_mpc);
    R_x = blkdiag(R_x,Zono.Rx_mpc);
end

R_u = blkdiag(R_u,Zono.Ru);
R_x = blkdiag(R_x,Zono.Rterm);

nsu = size(Zono.Ru,2)*N;
nsx = size(R_x,2);

S2Sx = [zeros(nsx,nsu),eye(nsx)];
S2Su = [eye(nsu),zeros(nsu,nsx)];

% S2Sx = [zeros(nsx,nsu),eye(nsx),zeros(nsx,n)];
% S2Su = [eye(nsu),zeros(nsu,nsx),zeros(nsu,n)];
% S2Sz = [zeros(n,nsu),zeros(n,nsx),eye(n)];
    
    
if ~isfield(prev_Data,'Linv')
    %% Nominal model
    
    prev_Data.Xnominal = Xk_1;
    
    
    %% define cost
    R_hat = [];
    Q_hat = [];
    dR_hat = [];
    dRmove = eye(N*nu)+diag(-1*ones((N-1)*nu,1),-nu);

    for i=1:N-1
        R_hat = blkdiag(R_hat,R);
        dR_hat = blkdiag(dR_hat,dR);
        Q_hat = blkdiag(Q_hat,C'*Q*C);
    end
    R_hat = blkdiag(R_hat,R);
    dR_hat = blkdiag(dR_hat,dR);
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
    
    H = S2Sx'*R_x'*Q_hat*R_x*S2Sx + S2Su'*R_u'*(R_hat + dR_hat)*R_u*S2Su;
    H = 0.00000001*eye(size(H,1)) + H;
    
    %% zono constraints
    Aineq = zeros((nsu+nsx)*2,(nsu+nsx));
    Bineq = ones((nsu+nsx)*2,1);

    for i=1:(nsu+nsx)
        Aineq(2*i-1,i)= 1;
        Aineq(2*i,i)= -1;
    end
    
    Aeq = R_x*S2Sx - Cu*R_u * S2Su ;
    
    

    
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
    
    %% modele nominal 
    
    prev_Data.Xnominal = A*prev_Data.Xnominal + Bu*prev_Data.u_prec + Bv*(Vm(1:nv)+ prev_Data.Vm);
    
    
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


Beq = - c_x + M*(prev_Data.Xnominal-x0) + Cv*Vm_pred + Cu*c_u;




U2u0 = zeros(nu,nu*N);
U2u0(:,1:nu) = eye(nu);

F = S2Sx'*R_x'*Q_hat*c_x + S2Su'*R_u'*(R_hat+dR_hat)*c_u - S2Su'*R_u'*U2u0'*dR*u_prec;


opt = mpcActiveSetOptions;
opt.UseHessianAsInput = false;

%[U,~,iA] = quadprog(Linv,F,Ac,Bc);
[S,info,iA] = mpcActiveSetSolver(Linv,F,Aineq,Bineq,Aeq,Beq,iA,opt);
prev_Data.iA = iA;

if info==-1 || info==-2
    disp("problem")
end

U = zeros(nu,N);
for i=1:N
    U(:,i) = Zono.Ru_mpc*S(nsu/N*(i-1)+1:nsu/N*i) + Zono.cu_mpc;
end



%% Intégrateur

% if ~isfield(prev_Data,'V')
%     prev_Data.V  = C*prev_Data.Xnominal - Ref;
% else
%     prev_Data.V = prev_Data.V + (C*prev_Data.Xnominal - Ref);
% end

prev_Data.u_prec = U(:,1);% - Ki*prev_Data.V;

%% robuste
U(:,1) = prev_Data.u_prec - K*(Xk_1 - prev_Data.Xnominal);

Zu = zonotope([Zono.cu, Zono.Ru]);
if ~in(Zu,U(:,1))
    disp("oups")
end
%%
Xplot = zeros(n,N+1);
Xplot (:,1) = prev_Data.Xnominal;
Sx = S2Sx*S;
for i=2:N
    %Xplot (:,i) = Zono.Rx_mpc*Sx(12*i-23:12*i-12) + Zono.cx_mpc;
    Xplot (:,i) = Zono.Rx_mpc*Sx(7*i-13:7*i-7) + Zono.cx_mpc;
    %Xplot (:,i) = Zono.Rx*Sx(2*i-3:2*i-2) + Zono.cx;
end  
Xplot (:,N+1) = Zono.Rterm*Sx(64:end) + Zono.cterm+x0; 
plot(Xplot(1,:),Xplot(2,:));
hold on
plot(Xplot(1,1),Xplot(2,1),'ro');
plot(zonotope([Zono.cterm+x0,Zono.Rterm]));
plot(x0(1,1),x0(2,1),'gx');
grid on



% q1=0.2; %discrétisation de la commande des volets
% u(1:2)= q1 * round(u(1:2)/q1);
% q2=5; %discrétisation de la commande des radiateurs
%u(3:5)= q2 * round(u(3:5)/q2);

new_data = prev_Data;
end