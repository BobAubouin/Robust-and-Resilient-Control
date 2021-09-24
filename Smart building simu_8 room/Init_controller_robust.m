load("reduced_model_coloc",'Tt','drsys_tc_12') % data from reduc model

%Définition du système de prédiction LTI
MPCobj.Plant = drsys_tc_12; 

A = drsys_tc_12.A;
Bs = drsys_tc_12.B;
C = drsys_tc_12.C;


MPCobj.iMV = 3:size_input; %index maniulated variable
MPCobj.iUD = size_input+1:size_input+(size_disturbances-size_predicted_disturbances); %index unmeasured disturbance
MPCobj.iMD = size_input+(size_disturbances-size_predicted_disturbances) + 1:size_input+size_disturbances; %index measured disturbance
MPCobj.iMO = 1:size_outputs; % index output

Bu = Bs(:,MPCobj.iMV);
Bd = Bs(:,[MPCobj.iUD MPCobj.iMD]);
Bv = Bs(:,[MPCobj.iUD MPCobj.iMD]);

Pert = VD.Data - Vm.Data;
pert = reshape(Pert,[size(Pert,1),size(Pert,3)]);
Zw.c = (max(pert,[],2) + min(pert,[],2))/2;

Zw.R = diag((max(pert,[],2) - min(pert,[],2))/2);
Zd1 = zonotope([Zw.c  Zw.R ]);
%gras
Zw.c = [Zw.c; zeros(size(A,1),1)];
Zw.R = blkdiag(Zw.R,eye(2)*0.00001,eye(size(A,1)-2)*0.0001);
Bd2 = [Bd, eye(size(A,1))];
% load("op_point","Xcomplet","E")
% Ze1 = zonotope.enclosePoints(Xcomplet(1:50,:));
% Ze2 = zonotope.enclosePoints(Xcomplet(51:100,:));
% Ze3 = zonotope.enclosePoints(Xcomplet(101:end,:));
% Zw.c = [Zw.c; Ze1.center; Ze2.center; Ze3.center];
% Zw.R = blkdiag(Zw.R,Ze1.generators,Ze2.generators,Ze3.generators);
% Zd2 = zonotope([Zw.c  Zw.R ]);
% Zw.c = [Zw.c; mean(Xcomplet,2)];
% Zw.R = blkdiag(Zw.R,diag(max(Xcomplet,[],2)-min(Xcomplet,[],2)));
% Bd2 = [Bd, E];
Zd = zonotope([Zw.c  Zw.R ]);


%Zd2 = reduce(Bd2*Zd,'scott',2);
Zd2 = Bd2*Zd;

% Horizon de prediction
MPCobj.PredictionHorizon = 10; 

% poids d'optimization
MPCobj.R = diag([5.1 1.8 9.4 8.9 11.2 9.9 24.6])/250; %poids commande (eye(7)/100;%)
MPCobj.dR = 0*eye(size_input-2); %poids variation de commande
MPCobj.Q = eye(size_outputs); %poids sur les sorties

MPCobj.Q(4,4) = 0.0001;

% contrainte commande sous forme de zonotope
zono.cu = [75;75;75;75;75;50;75];
zono.Ru = diag([75;75;75;75;75;50;75]);
ZU = zonotope([zono.cu zono.Ru]);

% contrainte sorties sous forme de zonotope
Tmax = -2500;
Tmin = 2000;

zono.cx = ones(12,1)*mean([Tmax,Tmin]);
zono.Rx = diag(repmat(mean([Tmax,Tmin])-Tmin,12,1));
ZX = zonotope([zono.cx zono.Rx]);

%% Robustesse
% p =    [0.6785;
%         0.7556;
%         0.8201;
%         0.99;
%         0.99;
%         0.999;
%         0.999;
%         0.999;
%         0.5;
%         0.5;
%         0.5;    
%         0.5];
% MPCobj.K = place(A,Bu(:,3:end),p);
Qlqg = eye(12);
Qlqg(10,10) = 0.1;
%R = eye(7)/10;
MPCobj.K = dlqr(A,Bu,C'*MPCobj.Q*C+eye(12)*10,MPCobj.R); % Calcul du retour d'état pour un fonctionnement en dual mode

%Calcul d'un d'un ensemble RPI pour borner la différence entre le modéle réel et nominal
alpha = 0.01;
Zset = getRPIset(A-Bu*MPCobj.K, Bd2 , Zd , alpha);
%Zset = deleteAligned(Zset);
% Plot_Zono(Zset)
% Plot_Zono((A-Bu*MPCobj.K)*Zset+Bd2*Zd,[],'r')
Zset = reduce(Zset,'combastel',10);
zono.cz = Zset.Z(:,1);
zono.Rz = Zset.Z(:,2:end);

% Calcul des contrainte restreinte pour le MPC
if length(zono.cu)>1
    Zumpc = ZU-(-MPCobj.K*Zset);
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
RHS = eye(12) + MPCobj.K'*MPCobj.R*MPCobj.K + (eye(size(A,2))-Phi)'*MPCobj.K'*MPCobj.dR*MPCobj.K*(eye(size(A,2))-Phi);
MPCobj.P = dlyap(Phi',RHS);

%% Initial condition and final objective
Mat = inv(eye(12)-A)*Bu;

Zref = Mat*ZU_mpc + Bv*Vmd2(:,1);
Ref = ones(8,1)*20;

% plot(Zref)
% hold on
% plot(Ref(1),Ref(2),'*')

if ~any(Ref)
    x0 = zeros(n,1);
    u0 = zeros(nu,1);
else
    cvx_begin
        variables x0(12,1) Su0(size(zono.Ru_mpc,2),1) r(size_outputs,1)
        minimize((r-Ref)'*MPCobj.Q*(r-Ref) + (zono.Ru_mpc*Su0 + zono.cu_mpc)'*MPCobj.R*(zono.Ru_mpc*Su0 + zono.cu_mpc))
        subject to
            [A-eye(12) Bu Bv; C zeros(size_outputs,size_input-2) zeros(size_outputs,6)] * [x0;zono.Ru_mpc*Su0 + zono.cu_mpc;Vmd2(:,1)] == [zeros(12,1);r];
            Su0 <= ones(size(zono.Ru_mpc,2),1);
            Su0 >= -ones(size(zono.Ru_mpc,2),1);
    cvx_end
end
u0 = zono.Ru_mpc*Su0 + zono.cu_mpc;
int.u0 = u0;
int.x0 = x0;
int.vm = Vmd2(:,1);

disp("Faisble ?")
disp(in(ZU_mpc,u0))
%% Set final

[~,~,Termset] = TerminalSet(MPCobj.K,MPCobj.P,ZU_mpc-u0,ZX_mpc-int.x0);
zono.cterm = Termset.Z(:,1);
zono.Rterm = Termset.Z(:,2:end);
X10 = InverseReach_inner(Termset,A,Bu,ZU_mpc-u0,ZX_mpc-x0,MPCobj.PredictionHorizon,0)+x0;
%on le met dans l'objet MPC
MPCobj.zono = zono;
% cvx_begin
%     variables s0(size(X10.Z,2)-1,1)
%     minimize((X10.Z(:,1)+X10.Z(:,2:end)*s0)'*eye(12)*(X10.Z(:,1)+X10.Z(:,2:end)*s0) + 0.000001*s0'*eye(size(X10.Z,2)-1)*s0)
%     subject to
%         s0 <= ones(size(X10.Z,2)-1,1);
%         s0 >= - ones(size(X10.Z,2)-1,1);
% cvx_end
% 
% x0 = X10.Z(:,1)+X10.Z(:,2:end)*s0/1.5;
disp("Initialisation terminée !")
%% plot
% dim = [2 3];
% % calcul de la zone d'attraction
% figure()
% plot(X10,dim,'FaceColor',[0 .6 0],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
% hold on
% 
% grid on
% %axis equal
% 
% 
% plot(Zset + int.x0 ,dim,'FaceColor',[0 0 0.8],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
% plot(Termset + int.x0 ,dim,'FaceColor',[0.6 .1 0],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
% plot(x0(dim(1)),x0(dim(2)),'*','LineWidth',2)
% plot(int.x0(1),int.x0(2),'*','LineWidth',2)
% legend("Initial feasible set",'$Z$',"$\hat{X}_f$","Starting point","Final point",'Interpreter','latex');
% ylabel('$x_2$','Interpreter','latex');
% xlabel('$x_1$','Interpreter','latex');
% %plot(ZX)
