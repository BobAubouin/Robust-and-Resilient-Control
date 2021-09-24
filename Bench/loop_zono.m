
init_scenar_ext;
init_plant;
%Init_full_controller;
Init_controller;
init_attacker_policy;
load('full_li_model_coloc.mat')
addpath("./test_simu")
%% init
tic
choix_scenario = 6; %1 = pas d'attaque, 2 = DoS, 3 = biais, 4 = upper sat, 5 = replay, 6 false prediction

ordre = size(MPCobj.Model.Plant.A,2);

ref = 22*ones(8,1);
T = 0:Ts:259200;
N = length(T);
U = zeros(size_input,N);
dU = zeros(size_input,N);
Y = zeros(size_outputs,N);
dY = zeros(size_outputs,N);
Y(:,1) = x0(1:size_outputs);
dY(:,1) = x0(1:size_outputs);
X = zeros(size_state,N);
X(:,1) = x0;
Xo = zeros(ordre,N);
Xo(:,1) = Tt(1:ordre,:)*x0;



BIY = variable0.BIY;
BIU = variable0.BIU;
Replay_data = variable0.RAYsave;

U_save = [];
mpc_data = struct;

for i = 1:N-1
   t = T(i);
   [X_k1,y] = Bilin(X(:,i),dU(:,i),VD.Data(:,:,i),A_k,Bu_k,Bv_k,Bxu_k,Bvu_k,C_k);
   %X_k1 = syslin.A*X(:,i) + syslin.B*[dU(:,i);VD.Data(:,:,i)];
   %y = syslin.C*X_k1;
   X(:,i+1) = X_k1;
   Y(:,i+1) = y;
   
   [dy,BIY,Replay_data,pred] = attacker_Y(y,Vm.Data(:,:,i:end),t,BIY,Replay_data,choix_scenario,Data,Ts);
   dY(:,i+1) = dy;
   obs = Tt(1:ordre,:)*X_k1;
   [u_mais,mpc_data,xo] = MPC_zono(MPCobj,dy,obs,ref,pred,U(:,i),L,mpc_data,zono);
   u = u_mais;
   Xo(:,i+1) = xo;
   U(:,i+1) = u;
   
   [du,BIU] = attacker_U(U(:,i+1),dU(:,i),t,BIU,choix_scenario,Data,Ts);
    dU(:,i+1) = du;
end
cost_bob = cost_calcul(Y,U,MPCobj);
toc
%% plot
errors = 0;
outputs = 1;
controls = 1;
observer = 0;
plot_simu_script;

%% save last run:

Y2 = Y;
dU2 = dU;
Xo2 = Xo;
disp("saved")

%%
errors = 0;
outputs = 1;
controls = 1;
observer = 0;
plot_compare_last_run

%%
figure()
p = [cost_matlab;cost_bob];
bar(p,'stacked')
legend("Pursuit error","Command power","Delta Command")
title("Cost of the trajectory")
xlabel("Matlab      vs      Bob")
ylabel("Total cost")
