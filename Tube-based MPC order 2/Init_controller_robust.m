%Définition du système de prédiction LTI
MPCobj.Plant = ss(A,[B,D],C,[],1); 

MPCobj.iMV = 1:size_input; %index maniulated variable
MPCobj.iUD = size_input+1:size_input+(size_disturbances-size_predicted_disturbances); %index unmeasured disturbance
MPCobj.iMD = size_input+(size_disturbances-size_predicted_disturbances) + 1:size_input+size_disturbances; %index measured disturbance
MPCobj.iMO = 1:size_outputs; % index output

% Horizon de prediction
MPCobj.PredictionHorizon = 10; 

% poids d'optimization
MPCobj.R = 0.1; %poids commande
MPCobj.dR = 0; %poids variation de commande
MPCobj.Q = eye(2); %poids sur les sorties

% contrainte commande sous forme de zonotope
zono.cu = 0;
zono.Ru = diag(1);
ZU = zonotope([zono.cu zono.Ru]);

% contrainte sorties sous forme de zonotope
Tmax = 5;
Tmin = -5;

zono.cx = ones(2,1)*mean([Tmax,Tmin]);
zono.Rx = diag(repmat(mean([Tmax,Tmin])-Tmin,2,1));
ZX = zonotope([zono.cx zono.Rx]);

%% Robustesse

MPCobj.K = dlqr(A,B,MPCobj.Q,MPCobj.R); % Calcul du retour d'état pour un fonctionnement en dual mode
%Calcul d'un d'un ensemble RPI pour borner la différence entre le modéle réel et nominal
alpha = 0.01;
Zset = getRPIset(A-B*MPCobj.K, D , Zd , alpha);
Zset = deleteAligned(Zset);
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

Zxmpc = ZX-Zset;
zono.cx_mpc = Zxmpc.Z(:,1);
zono.Rx_mpc = Zxmpc.Z(:,2:end);


%% poids final pour simulé un horizon de prédiction infini

Phi = A-B*MPCobj.K; 
RHS = C'*MPCobj.Q*C + MPCobj.K'*MPCobj.R*MPCobj.K + (eye(size(A,2))-Phi)'*MPCobj.K'*MPCobj.dR*MPCobj.K*(eye(size(A,2))-Phi);
MPCobj.P = dlyap(Phi',RHS);

%% Set final

ZU_mpc = zonotope([zono.cu_mpc zono.Ru_mpc]);
ZX_mpc = zonotope([zono.cx_mpc zono.Rx_mpc]);
[~,~,Termset] = TerminalSet(MPCobj.K,MPCobj.P,ZU_mpc,ZX_mpc);
zono.cterm = Termset.Z(:,1);
zono.Rterm = Termset.Z(:,2:end);



% calcul de la zone d'attraction
X10 = InverseReach(Termset,A,B,ZU,MPCobj.PredictionHorizon,0);
plot(X10,[1 2],'FaceColor',[0 .6 0],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
hold on
grid on
axis equal
plot(x0(1),x0(2),'*','LineWidth',2);
plot(Termset,[1 2],'FaceColor',[0.6 .1 0],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
plot(Zset,[1 2],'FaceColor',[0 0 0.8],'Filled',true,'FaceAlpha',1,'EdgeColor','none');
legend("Initial feasible set","Starting point","$\hat{X}_f$",'$Z$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
xlabel('$x_1$','Interpreter','latex')
%plot(ZX)
%% Zono atteignable

Mat = inv(eye(2)-A)*B;
zono.c0 = Mat*zono.cu;
zono.R0 = Mat*zono.Ru;
PolyU = mptPolytope(ZU_mpc);
PolyX = mptPolytope(ZX_mpc);
PolyTerm = mptPolytope(Termset);
MPCobj.Au = PolyU.P.A;
MPCobj.bu = PolyU.P.b;
MPCobj.Ax = PolyX.P.A;
MPCobj.bx = PolyX.P.b;
MPCobj.Aterm = PolyTerm.P.A;
MPCobj.bterm = PolyTerm.P.b;
%on le met dans l'objet MPC
MPCobj.zono = zono;