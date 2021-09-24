%% Parameters 
addpath(genpath('..'))
%systeme
A = [1   -1.2; 1.2 1.1];
B = [1 ; 0.5]; %NOK


P = [1.9385, 1.7088 ; 1.7088, 5.0552];
K = [-1.5718,- 0.2611];
D = [0.2 0.1; 0.4 0.1];
Zw.c = [0.02; 0];
Zw.R = [0.1 0; 0 0.005];
Zd = zonotope([Zw.c Zw.R]);
s = [1;1];

N= 10;
M =2;

Angle = deg2rad(0:1:180);
Forme = 0.1:0.1:1;
Forme = [Forme(1:end-1) flip(1./Forme)];
nb_angle = length(Angle);
Alpha = zeros(nb_angle,length(Forme));
Beta = zeros(nb_angle,length(Forme));
OK = zeros(nb_angle,length(Forme));

for l=1:length(Forme)
    H = [Forme(l) 0; 0 1/Forme(l)];
    for k=1:nb_angle
        H1 = rot(Angle(k))'*H*rot(Angle(k));
        [alpha,ok] = get_alpha(A,B,D,K,N,M,Zw,H1);
        OK(k,l) = ok;
        Alpha(k,l) = alpha;        
    end
end

% mesh(rad2deg(Angle),log(Forme),log(Alpha.*OK)')
% xlabel("orientation angle \theta (°)")
% ylabel("shape factor log(l)")
% zlabel("log(\alpha)")

%%
[alph_opt,i] = min(Alpha,[],'all','linear');
[angle_opt_i,Forme_opt_i] = ind2sub([nb_angle length(Forme)],i);
angle_opt = Angle(angle_opt_i);
H1_opt = rot(angle_opt)'*[Forme(Forme_opt_i) 0; 0 1/Forme(Forme_opt_i)]*rot(angle_opt);
ZA_opt = zonotope(ellipsoid(inv(H1_opt)*alph_opt),5,'i:norm');

for l=1:length(Forme)
    H = [Forme(l) 0; 0 1/Forme(l)];
    for k=1:nb_angle
        H2 = rot(Angle(k))'*H*rot(Angle(k));   
        Beta(k,l) = get_beta(A,B,D,K,N,M,Zw,H2,alph_opt,H1_opt);
    end
end
% mesh(rad2deg(Angle),log(Forme),log(Beta)')
% xlabel("orientation angle \theta (°)")
% ylabel("shape factor log(l)")
% zlabel("log(\alpha)")


[beta_opt,i] = min(Beta,[],'all','linear');
[angle_opt_i,Forme_opt_i] = ind2sub([nb_angle length(Forme)],i);
angle_opt = Angle(angle_opt_i);
H2_opt = rot(angle_opt)'*[Forme(Forme_opt_i) 0; 0 1/Forme(Forme_opt_i)]*rot(angle_opt);

ZB_opt = zonotope(ellipsoid(inv(H2_opt)*beta_opt),5,'i:norm');

% plot(ZA_opt,dim,'FaceColor',[0 0 0.8],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
% hold on
% plot(ZB_opt,dim,'FaceColor',[0 0 0.8],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
%%
MPCobj.Plant.A = A;
MPCobj.Plant.B = [B D];
MPCobj.Plant.C = eye(2);


MPCobj.iMV = 1; %index maniulated variable
MPCobj.iUD = 2:3; %index unmeasured disturbance
MPCobj.iMD = []; %index measured disturbance
MPCobj.iMO = 1:2; % index output
% Horizon de prediction
MPCobj.PredictionHorizon = 10; 

% poids d'optimization
MPCobj.R = 0.1; %poids commande
MPCobj.dR = 0; %poids variation de commande
MPCobj.Q = eye(2); %poids sur les sorties
MPCobj.P = P;
MPCobj.K = -K;
% contrainte commande sous forme de zonotope
zono.cu = [0];
zono.Ru = diag(110);
ZU = zonotope([zono.cu zono.Ru]);

% contrainte sorties sous forme de zonotope

zono.cx = zeros(2,1);
zono.Rx = diag([80 80]);
ZX = zonotope([zono.cx zono.Rx]);

%% calucl des ocntraintes restreintes
if length(zono.cu)>1
    Zumpc = ZU-(-MPCobj.K*ZA_opt);
    Zumpc2 = ZU-(-MPCobj.K*ZB_opt);
    zono.cu_mpc = Zumpc.Z(:,1);
    zono.Ru_mpc = Zumpc.Z(:,2:end);
    
else
    zono.cu_mpc = zono.cu + MPCobj.K*ZA_opt.Z(:,1);
    zono.Ru_mpc = abs(zono.Ru)-abs(MPCobj.K)*sum(abs(ZA_opt.Z(:,2:end)),2);
    zono.cu_mpc2 = zono.cu + MPCobj.K*ZB_opt.Z(:,1);
    zono.Ru_mpc2 = max(0,abs(zono.Ru)-abs(MPCobj.K)*sum(abs(ZB_opt.Z(:,2:end)),2));
end

Zxmpc = ZX-ZB_opt;
zono.cx_mpc = Zxmpc.Z(:,1);
zono.Rx_mpc = Zxmpc.Z(:,2:end);
ZU_mpc = zonotope([zono.cu_mpc zono.Ru_mpc]);
ZU_mpc2 = zonotope([zono.cu_mpc2 zono.Ru_mpc2]);
ZX_mpc = zonotope([zono.cx_mpc zono.Rx_mpc]);


%% Terminal and initial set

[~,~,Termset] = TerminalSet(K,MPCobj.P,ZU_mpc2,ZX_mpc);

zono.cterm = Termset.Z(:,1);
zono.Rterm = Termset.Z(:,2:end);
           
Xn_0 = InverseReach_1(Termset,A,B,ZU_mpc2,N,0);
Xn_1 = InverseReach_1(Termset,A,B,ZU_mpc2,N,1);
Xn_2 = InverseReach_1(Termset,A,B,ZU_mpc2,N,2);
x0 = randPoint(Xn_2);
MPCobj.zono = zono;
%% plot
dim = [1 2];
figure()
plot(Xn_0,dim,'FaceColor',[0 .6 0],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
hold on
plot(Xn_1,dim,'FaceColor',[0 .6 0.4],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
plot(Xn_2,dim,'FaceColor',[0 .6 0.8],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
grid on
%axis equal


plot(ZA_opt,dim,'FaceColor',[0 0 0.8],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
plot(Termset,dim,'FaceColor',[0.6 .1 0],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
% plot(x0(dim(1)),x0(dim(2)),'*','LineWidth',2)
% plot(0,0,'*','LineWidth',2)
legend("Initial feasible set (without DoS)","Initial feasible set(1 DoS) ","Initial feasible set(2 DoS)",'$\mu$RPI',"terminal set",'Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
xlabel('$x_1$','Interpreter','latex');
title("Sets invloved in the control scheme")
%% compare N-resilient, trajecory bound and RPI
figure()

plot(ZB_opt,dim,'FaceColor',[0.9 0.3 0.5],'Filled',true,'FaceAlpha',1,'EdgeColor','none');
hold on
plot(ZA_opt,dim,'FaceColor',[0 0.4 0.2],'Filled',true,'FaceAlpha',1,'EdgeColor','none');
alpha = 0.01;
RPI = getRPIset(A+B*K, D , Zd , alpha);
RPI = reduce(RPI,'combastel',10);
plot(10*RPI,[1 2],'FaceColor',[0 0 0.6],'Filled',true,'FaceAlpha',1,'EdgeColor','none');
legend("Trajectory bound (Z^+)","N-resilient set (Z)","RPI set (P)");
title("Compare N-resilient, trajecory bound and RPI")
%% test that the trajectory in Z will enter in the RPI
figure()

p1 = plot(ZA_opt,[1 2],'r');
hold on
Z = ZA_opt;
for i=1:7
    Z = (A+B*K)*Z + D*Zd;
    p2 = plot(Z);
end

disp("enter in the RPI ?")
disp(in(RPI,Z))
p3 = plot(RPI,[1 2],'g');
h = [p1 p2 p3];
legend(h,"N-Resilient set","Reachable set without DoS","RPI set")
title("Test that the trajectory in Z will enter in the RPI")
xlabel("z_1")
ylabel("z_2")