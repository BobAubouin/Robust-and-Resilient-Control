%% Parameters 
addpath(genpath('..'))
%systeme
A = [1.05   0.5; -0.6 1.1];
B = [1 ; 1]; %NOK

R = 0.1;
Q = eye(2);

K = -dlqr(A,B,Q,R);
Phi = A+B*K; 
RHS = Q + K'*R*K;
P = dlyap(Phi',RHS);


D = [1;0.5];
Zw.c = 0;
Zw.R = 0.1;
Zd = zonotope([Zw.c Zw.R]);
s = [1;1];

N= 10;

%%
alpha = 0.01;
RPI = getRPIset(A+B*K, D , Zd , alpha);
RPI = reduce(RPI,'combastel',10)*3;
%%
M =4;

Npoint = 100;
H = randn(2,Npoint);
for i=1:Npoint
    H(:,i) = H(:,i)/norm(H(:,i));
end
box = [1,0,-1,0;
       0,1,0,-1];
H=[H,box];
Bp = zeros(size(H,2),1);
tic
Bp(end+1) = 0;
[test_ok,Bnew] = TestMRPI(M,N,A,A+B*K,D*Zd,RPI,H,Bp);
toc
disp(['Nb itération : ',num2str(Bnew(end))])
Bnew = Bnew(1:end-1);

disp([' Le RPI est un N-RPI : ',num2str(test_ok)])
if test_ok
    Zplus = Polyhedron(H',Bnew);
    plot(mptPolytope(Zplus),[1,2],'r');
    hold on
    plot(RPI,[1,2],'b');
end
toc
ZA_opt = mptPolytope(RPI);
ZB_opt = mptPolytope(Zplus);

%%
MPCobj.Plant.A = A;
MPCobj.Plant.B = [B D];
MPCobj.Plant.C = eye(2);


MPCobj.iMV = 1; %index maniulated variable
MPCobj.iUD = 2; %index unmeasured disturbance
MPCobj.iMD = []; %index measured disturbance
MPCobj.iMO = 1:2; % index output
% Horizon de prediction
MPCobj.PredictionHorizon = 10; 

% poids d'optimization
MPCobj.R = R; %poids commande
MPCobj.dR = 0; %poids variation de commande
MPCobj.Q = Q; %poids sur les sorties
MPCobj.P = P;
MPCobj.K = -K;
% contrainte commande sous forme de zonotope
zono.cu = [0];
zono.Ru = diag(7);
ZU = zonotope([zono.cu zono.Ru]);

% contrainte sorties sous forme de zonotope

zono.cx = zeros(2,1);
zono.Rx = diag([20 20]);
ZX = zonotope([zono.cx zono.Rx]);

PX = mptPolytope(ZX);
PU = mptPolytope(ZU);
%% calucl des ocntraintes restreintes


V = vertices(ZB_opt);
P = mptPolytope(Polyhedron('V',V'));

Pu_mpc = PU-(-MPCobj.K*P);
Px_mpc = PX-ZB_opt;

V = vertices(Pu_mpc);
Pu_mpc = mptPolytope(Polyhedron('V',V','A',[1;-1],'b',abs(V(1))*ones(2,1)));
Pu_mpc = Pu_mpc.removeRedundancies;
Px_mpc = Px_mpc.removeRedundancies;

%% Terminal and initial set

[alpha,~,Termset] = TerminalSet_P(K,MPCobj.P,Pu_mpc,Px_mpc);
disp(['alpha = ',num2str(alpha)])
Pterm = mptPolytope(Termset);
           
Xn_0 = InverseReach_1(Termset,A,B,Pu_mpc,N,0);
Xn_1 = InverseReach_1(Termset,A,B,Pu_mpc,N,1);
Xn_4 = InverseReach_1(Termset,A,B,Pu_mpc,N,4);
%x0 = randPoint(Xn_2);
MPCobj.Au = Pu_mpc.P.A;
MPCobj.bu = Pu_mpc.P.b;
MPCobj.Ax = Px_mpc.P.A;
MPCobj.bx = Px_mpc.P.b;
MPCobj.Aterm = Pterm.P.A;
MPCobj.bterm = Pterm.P.b;
%% plot
dim = [1 2];
figure()
plot(Xn_0,dim,'FaceColor',[0 .6 0],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
hold on
plot(Xn_1,dim,'FaceColor',[0 .6 0.4],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
plot(Xn_4,dim,'FaceColor',[0 .6 0.8],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
grid on
%axis equal


plot(ZA_opt,dim,'FaceColor',[0 0 0.8],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
plot(Termset,dim,'FaceColor',[0.6 .1 0],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
plot(Px_mpc,dim);
% plot(x0(dim(1)),x0(dim(2)),'*','LineWidth',2)
% plot(0,0,'*','LineWidth',2)
legend("Initial feasible set (without DoS)","Initial feasible set(1 DoS) ","Initial feasible set(4 DoS)",'N-Resilient',"Terminal set",'Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
xlabel('$x_1$','Interpreter','latex');
%title("Sets invloved in the control scheme",'Interpreter','latex')
%% compare N-resilient, trajecory bound and RPI
figure()

plot(ZB_opt,dim,'FaceColor',[0.9 0.3 0.5],'Filled',true,'FaceAlpha',1,'EdgeColor','none');
hold on
plot(ZA_opt,dim,'FaceColor',[0 0.4 0.2],'Filled',true,'FaceAlpha',1,'EdgeColor','none');
alpha = 0.01;
RPI = getRPIset(A+B*K, D , Zd , alpha);
RPI = reduce(RPI,'combastel',10);
plot(RPI,[1 2],'FaceColor',[0 0 0.6],'Filled',true,'FaceAlpha',1,'EdgeColor','none');
legend("Trajectory bound (Z^+)","N-resilient set (Z)","RPI set (P)");
title("Compare N-RPI, trajecory bound and RPI")
