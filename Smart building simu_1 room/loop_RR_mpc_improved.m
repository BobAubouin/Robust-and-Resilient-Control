addpath(genpath(".."))
rng(1)
init_scenar_ext;
init_plant;
%Init_full_controller;
Init_controller_robust;
init_attacker_policy;

%% init
tic

choix_scenario = 1; %1 = pas d'attaque, 2 = DoS, 3 = biais, 4 = upper sat, 5 = replay, 6 false prediction

ordre = size(MPCobj.Plant.A,2);


Ts = 15*60; %discretization toute les 15 minutes
T = 0:Ts/3600:72; %temps en heure
N = length(T);
U = zeros(size_input,N);
U_mpc = zeros(size_input,N);
dU = zeros(size_input,N);
Y = zeros(size_outputs,N);
dY = zeros(size_outputs,N);
Y(:,1) = x0(1:size_outputs);
dY(:,1) = x0(1:size_outputs);
X = zeros(size_state,N);
X(:,1) = x0;
Xo = zeros(ordre,N);
Xo(:,1) = x0;

XN = zeros(ordre,N);
XN(:,1) = x0;

X_pred = zeros(ordre,N);
X_pred(:,1) = x0;

BIY = variable0.BIY;
BIU = variable0.BIU;
Replay_data = variable0.RAYsave;

U_save = [];


mpc_data = struct;
mpc_data.x0 = int.x0;
mpc_data.u0 = int.u0;
u_mais = [];

Ref = zeros(ordre,N);
Ref(:,1) = mpc_data.x0;
for i = 1:N-1
    if i==100
        mpc_data.x0 = int.x0_2;
        mpc_data.u0 = int.u0_2;
    end
   Ref(:,i+1)= mpc_data.x0;
   t = T(i);
   y = C*X(:,i);
   [dy,BIY,Replay_data,pred] = attacker_Y(y,Vm.Data(:,:,i:end),t,BIY,Replay_data,choix_scenario,Data,Ts); %Vm.Data(:,:,i:end)
   dY(:,i+1) = dy;
   obs = X(:,i);
   if Dos(i)
       if isempty(u_mais)
           u = 0;
       else
           nb_dos = nb_dos+1;
           %u = u_mais(nb_dos);
           %u = 0;
           %u = mpc_data.u0;
           u = U(:,i-1);
           disp(nb_dos)
           U_mpc(:,i) = u_mais(nb_dos);
           X_pred(:,i) = x_pred(:,nb_dos);
       end
   else
       nb_dos = 1;
       if i==1
           [u_mais,mpc_data,xn,x_pred] = MPC_ROBUSTE_improved(MPCobj,obs,Ref,pred,U(:,1),mpc_data);
       else
           [u_mais,mpc_data,xn,x_pred] = MPC_ROBUSTE_improved(MPCobj,obs,Ref,pred,U(:,i-1),mpc_data);
       end
       u = u_mais(:,1);
       U_mpc(:,i) = mpc_data.u_prec;
       X_pred(:,i) = x_pred(:,1);
   end
   
   Xo(:,i) = xn;
   XN(:,i) = xn;
   U(:,i) = u;
   
   
   dU(:,i) = u;
    
   %[X_k1,y] = Bilin(X(:,i),dU(:,i),VD.Data(:,:,i),A_k,Bu_k,Bv_k,Bxu_k,Bvu_k,C_k);
   X_k1 = A*X(:,i) + Bu*dU(:,i) + Bd*(VD.Data(:,:,i+1)-0*[0;Vm.Data(:,:,i+1)]);
   y = C*X_k1;
   X(:,i+1) = X_k1;
   Y(:,i+1) = y;
   disp(i)
end
toc

XN(:,N) = XN(:,N-1);
X_pred(:,N) = X_pred(:,N-1);
XNplus = X_pred + zono.cz + abs(zono.Rz)*ones(size(zono.Rz,2),N);
XNmoins = X_pred + zono.cz - abs(zono.Rz)*ones(size(zono.Rz,2),N);
cost_improved = cost_calcul(Y,U,MPCobj,mpc_data);
%% plot

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
