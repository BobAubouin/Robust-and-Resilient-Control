addpath(genpath(".."))
rng(1)
init_scenar_ext;
init_plant;
%Init_full_controller;
Init_controller_robust;
T = 0:Ts/3600:72; %temps en heure
init_attacker_policy;

%% init
tic

choix_scenario = 1; %1 = pas d'attaque, 2 = DoS, 3 = biais, 4 = upper sat, 5 = replay, 6 false prediction

ordre = size(MPCobj.Plant.A,2);

x0 = int.x0;
Ts = 15*60; %discretization toute les 15 minutes
N = length(T);
size_input = 7;
U = zeros(size_input,N);
U_mpc = zeros(size_input,N);
dU = zeros(size_input,N);
Y = zeros(size_outputs,N);
dY = zeros(size_outputs,N);
Y(:,1) = C*x0;
dY(:,1) = C*x0;
X = zeros(12,N);
X(:,1) = x0;
Xo = zeros(ordre,N);
XN = zeros(ordre,N);


BIY = variable0.BIY;
BIU = variable0.BIU;
Replay_data = variable0.RAYsave;

U_save = [];

mpc_data = struct;
mpc_data.x0 = int.x0;
mpc_data.u0 = int.u0;
for i = 1:N-1
   t = T(i);
   y = C*X(:,i);
   [dy,BIY,Replay_data,pred] = attacker_Y(y,Vm.Data(:,:,i:end) - Vm2.Data(:,:,i:end),t,BIY,Replay_data,choix_scenario,Data,Ts); %Vm.Data(:,:,i:end)
   dY(:,i+1) = dy;
   obs = X(:,i);
   if i==1
       [u_mais,mpc_data,xn] = MPC_ROBUSTE_improved(MPCobj,obs,Ref,pred,U(:,1),mpc_data,0);
   else
       [u_mais,mpc_data,xn] = MPC_ROBUSTE_improved(MPCobj,obs,Ref,pred,U(:,i-1),mpc_data,0);
   end
   
   u = u_mais(:,1);
   Xo(:,i) = xn;
   XN(:,i) = xn;
   U(:,i) = u;
   U_mpc(:,i) = mpc_data.u_prec;
   
   if i==1
       [du,BIU] = attacker_U(U(:,i),dU(:,i),t,BIU,choix_scenario,Data,Ts);
   else  
       [du,BIU] = attacker_U(U(:,i),dU(:,i-1),t,BIU,choix_scenario,Data,Ts);
   end
    dU(:,i) = du;
    
   %[X_k1,y] = Bilin(X(:,i),dU(:,i),VD.Data(:,:,i),A_k,Bu_k,Bv_k,Bxu_k,Bvu_k,C_k);
   X_k1 = A*X(:,i) + Bu*dU(:,i) + Bv*(VD.Data(:,:,i+1).*[1;1;1;1;1;1]);%-1*[0;0;0;Vm.Data(:,:,i+1)]);
   y = C*X_k1;
   X(:,i+1) = X_k1;
   Y(:,i+1) = y;
   disp(i)
end
toc

XN(:,N) = XN(:,N-1);
YNplus = C*XN + C*zono.cz + abs(C*zono.Rz)*ones(size(C*zono.Rz,2),N);
YNmoins = C*XN + C*zono.cz - abs(C*zono.Rz)*ones(size(C*zono.Rz,2),N);
cost_improved = cost_calcul(Y,U,MPCobj,mpc_data);
%% plot
Ref = mpc_data.x0;
dU(:,end)=dU(:,end-2);
% dU(:,end-1)=dU(:,end-2);
U_mpc(:,end)=U_mpc(:,end-2);
errors = 0;
outputs = 1;
controls = 1;
observer = 0;
zonoplot = 0;
scenar_plot = 0;
plot_simu_script;

%% save last run:

Y2 = Y;
dU2 = dU;
Xo2 = Xo;
XN2 = XN;
zono2 = zono;

disp("saved")

%%
errors = 0;
outputs = 0;
controls = 1;
observer = 0;
zonoplot = 1;
plot_compare_last_run

%%
figure()
p = [cost_robuste;cost_improved];
bar(p,'stacked')
legend("Pursuit error","Command power","Delta Command")
title("Cost of the trajectory")
xlabel("Standard      vs      Improved")
ylabel("Total cost")
