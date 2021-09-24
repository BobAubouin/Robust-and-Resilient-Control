
%Définition du système de prédiction LTI
C_k = eye(7);
size_outputs = 7;
MPCobj.Plant = ss(A_k,[Bu_k(:,2) , Bv_k],C_k,[],1); 

size_input = 1;
A = A_k;
C = C_k;


MPCobj.iMV = 1;%1:size_input; %index maniulated variable
MPCobj.iUD = size_input+1:size_input+(size_disturbances-size_predicted_disturbances); %index unmeasured disturbance
MPCobj.iMD = size_input+(size_disturbances-size_predicted_disturbances) + 1:size_input+size_disturbances; %index measured disturbance
MPCobj.iMO = 1:size_outputs; % index output

Bu = Bu_k(:,2);
Bd = Bv_k;
% 
Zw.c = [ 2.5; 0; 13];
Zw.R = diag([ 2.5; 1; 25]);
%gras
Zw.c = [Zw.c; zeros(size(A,1),1)];
Zw.R = blkdiag(Zw.R,eye(size(A,1))*0.0001);
Bd2 = [Bd, eye(size(A,1))];
Zd = zonotope([Zw.c  Zw.R ]);
% Pert = VD.Data - Vm.Data;
% pert = reshape(Pert,[size(Pert,1),size(Pert,3)]);
% Zw.c = (max(pert,[],2) + min(pert,[],2))/2;
% 
% Zw.R = diag((max(pert,[],2) - min(pert,[],2))/2);
% Zw.c = [Zw.c; zeros(size(A,1),1)];
% Zw.R = blkdiag(Zw.R,eye(2)*0.00001,eye(size(A,1)-2)*0.0001);
% Bd2 = [Bd, eye(size(A,1))];
% Zd = zonotope([Zw.c  Zw.R ]);


% Horizon de prediction
MPCobj.PredictionHorizon = 10; 

% poids d'optimization
MPCobj.R = 0.1*eye(1); %poids commande
MPCobj.dR = 0*eye(1); %poids variation de commande
MPCobj.Q = diag([100,0.1,0.1,0.1,0.1,0.1,0.1]); %poids sur les sorties

% contrainte commande sous forme de zonotope
zono.cu = [50];
zono.Ru = diag([50]);
ZU = zonotope([zono.cu zono.Ru]);

% contrainte sorties sous forme de zonotope
Tmax = 50;
Tmin = 0;

zono.cx = ones(size_state,1)*mean([Tmax,Tmin]);
zono.Rx = diag(repmat(mean([Tmax,Tmin])-Tmin,size_state,1));
ZX = zonotope([zono.cx zono.Rx]);

%% Robustesse

MPCobj.K = dlqr(A,Bu,MPCobj.Q/5,MPCobj.R); % Calcul du retour d'état pour un fonctionnement en dual mode
%Calcul d'un d'un ensemble RPI pour borner la différence entre le modéle réel et nominal
alpha = 0.01;
Zset = getRPIset(A-Bu*MPCobj.K, Bd2 , Zd , alpha);
%Zset = deleteAligned(Zset);
Zset = reduce(Zset,'combastel',5);
zono.cz = Zset.Z(:,1);
zono.Rz = Zset.Z(:,2:end);

% Calcul des contrainte restreinte pour le MPC
if length(zono.cu)>1
    Zumpc = ZU-MPCobj.K*Zset;
    zono.cu_mpc = Zumpc.Z(:,1);
    zono.Ru_mpc = Zumpc.Z(:,2:end);
else
    zono.cu_mpc = zono.cu + MPCobj.K*zono.cz;
    zono.Ru_mpc = abs(zono.Ru)-abs(MPCobj.K)*sum(abs(zono.Rz),2);
end

Zxmpc = ZX;
zono.cx_mpc = Zxmpc.Z(:,1);
zono.Rx_mpc = Zxmpc.Z(:,2:end);
ZU_mpc = zonotope([zono.cu_mpc zono.Ru_mpc]);
ZX_mpc = zonotope([zono.cx_mpc zono.Rx_mpc]);

%% poids final pour simulé un horizon de prédiction infini

Phi = A-Bu*MPCobj.K; 
RHS = MPCobj.Q+ MPCobj.K'*MPCobj.R*MPCobj.K + (eye(size_state)-Phi)'*MPCobj.K'*MPCobj.dR*MPCobj.K*(eye(size_state)-Phi);
MPCobj.P = dlyap(Phi',RHS);

%% Initial condition and final objective


Mat = inv(eye(size_state)-A)*Bu;

Zref = Mat*ZU_mpc + Bd(:,2:3)*Vmd2(:,1);
Ref = [20;18;13;ones(3,1)*15;18];
% plot(Zref)
% hold on
% plot(Ref(1),Ref(2),'*')

if ~any(Ref)
    x0 = zeros(n,1);
    u0 = zeros(nu,1);
else
    cvx_begin
        variables x0(size_state,1) u0(size_input,1) r(size_outputs,1)
        minimize((r-Ref)'*(MPCobj.Q)*(r-Ref))
        subject to
            [A-eye(size_state) Bu Bd(:,2:3); C zeros(size_outputs,size_input) zeros(size_outputs,size_disturbances-1)] * [x0;u0;Vmd2(:,1)] == [zeros(size_state,1);r];
    cvx_end
end
int.u0 = u0;
int.x0 = x0;
int.vm = Vmd2(:,1);

Ref2 = [18.5;18.5;14;ones(3,1)*16;18.5];
if ~any(Ref2)
    x0 = zeros(n,1);
    u0 = zeros(nu,1);
else
    cvx_begin
        variables x0(size_state,1) u0(size_input,1) r(size_outputs,1)
        minimize((r-Ref2)'*(MPCobj.Q*10)*(r-Ref2))
        subject to
            [A-eye(size_state) Bu Bd(:,2:3); C zeros(size_outputs,size_input) zeros(size_outputs,size_disturbances-1)] * [x0;u0;Vmd2(:,1)] == [zeros(size_state,1);r];
    cvx_end
end
int.u0_2 = u0;
int.x0_2 = x0;

%% Set final


[~,~,Termset] = TerminalSet(MPCobj.K,MPCobj.P,ZU_mpc-int.u0,ZX_mpc-int.x0) ;
zono.cterm = Termset.Z(:,1);
zono.Rterm = Termset.Z(:,2:end);
X10 = InverseReach_inner(Termset,A,Bu,ZU_mpc-int.u0,ZX_mpc-int.x0,MPCobj.PredictionHorizon,0)+int.x0;

cvx_begin
    variables s0(size(X10.Z,2)-1,1)
    minimize((X10.Z(:,1)+X10.Z(:,2:end)*s0)'*eye(size_state)*(X10.Z(:,1)+X10.Z(:,2:end)*s0) + 0.000001*s0'*eye(size(X10.Z,2)-1)*s0)
    subject to
        s0 <= ones(size(X10.Z,2)-1,1);
        s0 >= - ones(size(X10.Z,2)-1,1);
cvx_end


x0 = X10.Z(:,1)+X10.Z(:,2:end)*s0;
%x0 = randPoint(X10);

%% plot
dim = [1 2];
% calcul de la zone d'attraction
figure()
plot(X10,dim,'FaceColor',[0 .6 0],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
hold on

grid on
%axis equal


plot(Zset + int.x0 ,dim,'FaceColor',[0 0 0.8],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
plot(Termset + int.x0 ,dim,'FaceColor',[0.6 .1 0],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
plot(x0(dim(1)),x0(dim(2)),'*','LineWidth',2)
plot(int.x0(dim(1)),int.x0(dim(2)),'*','LineWidth',2)
legend("Initial feasible set",'$Z$',"$\hat{X}_f$","Starting point","Final point",'Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
xlabel('$x_1$','Interpreter','latex');
%plot(ZX)
%% Zono atteignable

%on le met dans l'objet MPC
tic
Termset = reduce(Termset,'combastel',3);
PolyU = mptPolytope(ZU_mpc);
PolyX = mptPolytope(ZX_mpc);
PolyTerm = mptPolytope(Termset);
MPCobj.Au = PolyU.P.A;
MPCobj.bu = PolyU.P.b;
MPCobj.Ax = PolyX.P.A;
MPCobj.bx = PolyX.P.b;
MPCobj.Aterm = PolyTerm.P.A;
MPCobj.bterm = PolyTerm.P.b;
MPCobj.zono = zono;
toc