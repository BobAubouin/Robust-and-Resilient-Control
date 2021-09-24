%% Parameters 
addpath(genpath('..'))
%systeme
load("discrete_model_solo") % data from BRCM

%full discrete systeme
A = B.building_model.discrete_time_model.A; % same for Bu,Bv,Bvu,Bxu
D = B.building_model.discrete_time_model.Bv;
B = B.building_model.discrete_time_model.Bu;
C = [1 0 0 0 0 0 0];

B = B(:,2); % on enlève les volets

Zw.c = [0];
Zw.R = diag([0.001]);
%gras
% Zw.c = [Zw.c; zeros(size(A,1),1)];
% Zw.R = blkdiag(Zw.R,eye(size(A,1))*0.0001);
% Bd2 = [D(:,2), eye(size(A,1))];
Zd = zonotope([Zw.c Zw.R]);

N= 10;

MPCobj.Plant.A = A;
MPCobj.Plant.B = [B D];
MPCobj.Plant.C = eye(7);


MPCobj.iMV = 1; %index maniulated variable
MPCobj.iUD = 2; %index unmeasured disturbance
MPCobj.iMD = 3:4; %index measured disturbance
MPCobj.iMO = 1; % index output
% Horizon de prediction
MPCobj.PredictionHorizon = N; 

% poids d'optimization
MPCobj.R = 0.1*eye(1); %poids commande
MPCobj.dR = 0*eye(1); %poids variation de commande
MPCobj.Q = diag([100,1,1,1,1,1,1]); %poids sur les sorties

% contrainte commande sous forme de zonotope
Pmax = 100;
zono.cu = [Pmax/2];
zono.Ru = diag([Pmax/2]);
ZU = zonotope([zono.cu zono.Ru]);

% contrainte sorties sous forme de zonotope
Tmax = 30;
Tmin = 10;

zono.cx = ones(7,1)*mean([Tmax,Tmin]);
zono.Rx = diag(repmat(mean([Tmax,Tmin])-Tmin,7,1));
ZX = zonotope([zono.cx zono.Rx]);
%% Calcul of the feedback

MPCobj.K = - dlqr(A,B,MPCobj.Q*10,MPCobj.R);
Phi = A+B*MPCobj.K; 
RHS = MPCobj.Q+ MPCobj.K'*MPCobj.R*MPCobj.K ;
MPCobj.P = dlyap(Phi',RHS); %terminal weight

%% Calcul of the invariant sets
% alpha = 0.01;
% Zset = getRPIset(A+B*MPCobj.K, Bd2 , Zd , alpha);
% Zset = reduce(Zset,'combastel',5);
% Approx_RPI = ellipsoid(Zset);
% 
% invH = Approx_RPI.Q;

%H = inv(invH);%MPCobj.P;%
H2 = MPCobj.P;
invH2 = inv(H2);
%invH = inv(H);
M = 2;
[alpha,ok] = get_alpha(A,B,D(:,2),MPCobj.K,N,M,Zw,H2);
disp('max Number attack ok:')
disp(ok)
EA = ellipsoid(invH2*alpha);
ZA_opt = zonotope(EA,100,'i:norm')*diag([2,2,2,4,2,2,2]);
Beta = get_beta(A,B,D(:,2),MPCobj.K,N,M,Zw,H2,alpha,H2);
EB = ellipsoid(invH2*Beta);
ZB_opt = zonotope(EB,100,'i:norm')*1.5;

dim = [1 2];
plot(ZA_opt,dim,'FaceColor',[0 0 0.8],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
hold on
plot(ZB_opt,dim,'FaceColor',[0 0 0.8],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');

%% calucl des ocntraintes restreintes
if length(zono.cu)>1
    Zumpc = ZU-(MPCobj.K*ZA_opt);
    Zumpc2 = ZU-(MPCobj.K*ZB_opt);
    zono.cu_mpc = Zumpc.Z(:,1);
    zono.Ru_mpc = Zumpc.Z(:,2:end);
    
else
    zono.cu_mpc = zono.cu - MPCobj.K*ZA_opt.Z(:,1);
    zono.Ru_mpc = abs(zono.Ru)-abs(MPCobj.K)*sum(abs(ZA_opt.Z(:,2:end)),2);
    zono.cu_mpc2 = zono.cu - MPCobj.K*ZB_opt.Z(:,1);
    zono.Ru_mpc2 = max(0,abs(zono.Ru)-abs(MPCobj.K)*sum(abs(ZB_opt.Z(:,2:end)),2));
end

Zxmpc = ZX;
zono.cx_mpc = Zxmpc.Z(:,1);
zono.Rx_mpc = Zxmpc.Z(:,2:end);
ZU_mpc = zonotope([zono.cu_mpc zono.Ru_mpc]);
ZU_mpc2 = zonotope([zono.cu_mpc2 zono.Ru_mpc2]);
ZX_mpc = zonotope([zono.cx_mpc zono.Rx_mpc]);



%%
C = eye(7);
Ref = [20.8036
   20.8100
   10.6292
   10.6292
   10.6292
   10.6292
   20.8131];
if ~any(Ref)
    x0 = zeros(n,1);
    u0 = zeros(nu,1);
else
    cvx_begin
        variables x0(7,1) u0(1,1) r(7,1)
        minimize((r-Ref)'*eye(7)*(r-Ref))
        subject to
            [A-eye(7) B D(:,2:3); C zeros(7,1) zeros(7,2)] * [x0;u0;Vmd2(:,1)] == [zeros(7,1);r];
    cvx_end
end
int.u0 = u0;
int.x0 = x0;
int.vm = Vmd2(:,1);
disp("faisable ?")
disp(in(ZU_mpc2,u0));
%% Terminal and initial set

[~,~,Termset] = TerminalSet(-MPCobj.K,MPCobj.P,ZU_mpc2-int.u0,ZX_mpc-int.x0);

zono.cterm = Termset.Z(:,1);
zono.Rterm = Termset.Z(:,2:end);
           
Xn_0 = InverseReach_1(Termset,A,B,ZU_mpc2,N,0);
Xn_1 = InverseReach_1(Termset,A,B,ZU_mpc2,N,1);
Xn_2 = InverseReach_1(Termset,A,B,ZU_mpc2,N,2);
x0 = randPoint(Xn_2);

MPCobj.K = - MPCobj.K;
MPCobj.zono = zono;
%% plot
% dim = [1 2];
% figure()
% plot(Xn_0,dim,'FaceColor',[0 .6 0],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
% hold on
% plot(Xn_1,dim,'FaceColor',[0 .6 0.4],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
% plot(Xn_2,dim,'FaceColor',[0 .6 0.8],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
% grid on
% %axis equal
% 
% 
% plot(ZA_opt,dim,'FaceColor',[0 0 0.8],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
% plot(Termset,dim,'FaceColor',[0.6 .1 0],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
% % plot(x0(dim(1)),x0(dim(2)),'*','LineWidth',2)
% % plot(0,0,'*','LineWidth',2)
% legend("Initial feasible set (without DoS)","Initial feasible set(1 DoS) ","Initial feasible set(2 DoS)",'$\mu$RPI',"terminal set","Starting point","Final point",'Interpreter','latex');
% ylabel('$x_2$','Interpreter','latex');
% xlabel('$x_1$','Interpreter','latex');

%% 
% 
% plot(ZA_opt,dim,'FaceColor',[0 0.3 0.8],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
% hold on
% plot(ZB_opt,dim,'FaceColor',[0.9 0 0.8],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
% alpha = 0.01;
% RPI = getRPIset(A+B*K, D , Zd , alpha);
% RPI = reduce(RPI,'combastel',10);
% plot(10*RPI,[1 2],'FaceColor',[0 0 0.6],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');

