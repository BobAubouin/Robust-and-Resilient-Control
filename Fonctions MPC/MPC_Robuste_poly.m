function [U,new_data,Xk] = MPC_Robuste_poly(mpc,Xk_1,Ref,Vm,u_prec,prev_Data,model)

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
        Bv = [B(:,mpc.iMD), Bd];
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
u0 = prev_Data.u0;
x0 = prev_Data.x0;


    
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
    
    H = Cu'*Q_hat*Cu +  R_hat + dR_hat;
    H = 0.00000001*eye(size(H,1)) + H;
   
    
    

    
    L = chol(H,'lower');
    Linv = L\eye(size(H,1));
    
    
    prev_Data.Linv = Linv;
    prev_Data.H = H;
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
    iA = prev_Data.iA;
    M = prev_Data.M;
    Cv = prev_Data.Cv;
    Cu = prev_Data.Cu;
    Q_hat = prev_Data.Q_hat;
    R_hat = prev_Data.R_hat;
    dR_hat = prev_Data.dR_hat;
end



%%
U2u0 = zeros(nu,nu*N);
U2u0(:,1:nu) = eye(nu);

Beta = M*(prev_Data.Xnominal-x0)+Cv*Vm_pred;

F = Cu'*Q_hat*Beta - U2u0'*dR*u_prec;

%% polytopic constraints
Au = mpc.Au;
bu = mpc.bu;
Ax = mpc.Ax;
bx = mpc.bx;
Aterm = mpc.Aterm;
bterm = mpc.bterm;

AU=[];
AX =[];
for i=1:N-1
    AU = blkdiag(AU,Au);
    AX = blkdiag(AX,Ax);
end
AU = blkdiag(AU,Au);
BU = repmat(bu-Au*u0,N,1);
AX = blkdiag(AX, Aterm);
AX2=AX*Cu;
BX = repmat(bx-Ax*x0,N-1,1);
BX = [BX;bterm];
BX = BX-AX*Beta;

Aineq = [AU;AX2];
Bineq = [BU;BX];

if ~isfield(prev_Data,'iA')
    iA = false(size(Bineq,1),1);
    prev_Data.iA = iA;
end

%% OPtimization
opt = mpcActiveSetOptions;
opt.UseHessianAsInput = false;

%[U,~,iA] = quadprog(Linv,F,Ac,Bc);
[S,info,iA] = mpcActiveSetSolver(Linv,F,Aineq,Bineq,zeros(0,N*nu),zeros(0,1),iA,opt);
prev_Data.iA = iA;

if info==-1 || info==-2
    disp("problem")
end
U=[];
for i=1:N
    U(:,i)=S(i*nu:(i+1)*nu-1)+u0;
end



%% Intégrateur

% if ~isfield(prev_Data,'V')
%     prev_Data.V  = C*prev_Data.Xnominal - Ref;
% else
%     prev_Data.V = prev_Data.V + (C*prev_Data.Xnominal - Ref);
% end

prev_Data.u_prec = U(1);% - Ki*prev_Data.V;

%% robuste
U(:,1) = prev_Data.u_prec - K*(Xk_1 - prev_Data.Xnominal);
% Zono = mpc.zono;
% Zu = zonotope([Zono.cu, Zono.Ru]);
% if ~in(Zu,U(:,1))
%     disp("oups")
% end
%%
% Xplot = zeros(n,N+1);
% Xplot (:,1) = prev_Data.Xnominal;
% Sx = Beta+Cu*S;
% for i=2:N+1 
%     %Xplot (:,i) = Zono.Rx_mpc*Sx(12*i-23:12*i-12) + Zono.cx_mpc;
%     %Xplot (:,i) = Sx(7*i-13:7*i-7)+x0;
%     Xplot (:,i) = Sx(2*i-3:2*i-2)+x0;
% end  
% plot(Xplot(1,:),Xplot(2,:))
% hold on
% plot(Xplot(1,1),Xplot(2,1),'ro')
% plot(x0(1,1),x0(2,1),'gx')
% grid on



% q1=0.2; %discrétisation de la commande des volets
% u(1:2)= q1 * round(u(1:2)/q1);
% q2=5; %discrétisation de la commande des radiateurs
%u(3:5)= q2 * round(u(3:5)/q2);

new_data = prev_Data;
end