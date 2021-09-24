addpath(genpath(".."))
init_plant;
%Init_full_controller;
Init_controller_robust;
T = 0:Ts:100;
init_attacker_policy;

%% init
tic
rng(1)
choix_scenario = 6; %1 = pas d'attaque, 2 = DoS, 3 = biais, 4 = upper sat, 5 = replay, 6 false prediction

ordre = size(MPCobj.Plant.A,2);

ref = zono.c0-zono.R0*[0.6];
ref = [0; 0];

T = 0:Ts:100;
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
Xo(:,1) = x0;

XN = zeros(ordre,N);
XN(:,1) = x0;

BIY = variable0.BIY;
BIU = variable0.BIU;
Replay_data = variable0.RAYsave;

U_save = [];
mpc_data = struct;

for i = 1:N-1
   t = T(i);
   y = C*X(:,i);
   [dy,BIY,Replay_data,pred] = attacker_Y(y,[],t,BIY,Replay_data,choix_scenario,Data,Ts);
   dY(:,i+1) = dy;
   obs = X(:,i);
   if i==1
       [u_mais,mpc_data,xn] = MPC_ROBUSTE_improved(MPCobj,obs,ref,pred,U(:,1),mpc_data);
   else
       [u_mais,mpc_data,xn] = MPC_ROBUSTE_improved(MPCobj,obs,ref,pred,U(:,i-1),mpc_data);
   end
   
   u = u_mais;
   Xo(:,i) = xn;
   XN(:,i) = xn;
   U(:,i) = u;
   
   if i==1
       [du,BIU] = attacker_U(U(:,i),dU(:,i),t,BIU,choix_scenario,Data,Ts);
   else  
       [du,BIU] = attacker_U(U(:,i),dU(:,i-1),t,BIU,choix_scenario,Data,Ts);
   end
    dU(:,i) = du;
    
   %[X_k1,y] = Bilin(X(:,i),dU(:,i),VD.Data(:,:,i),A_k,Bu_k,Bv_k,Bxu_k,Bvu_k,C_k);
   X_k1 = A*X(:,i) + B*dU(:,i) + D*(Zw.c + Zw.R*(rand(2,1)*2-1));
   y = C*X_k1;
   X(:,i+1) = X_k1;
   Y(:,i+1) = y;

   disp(i)
end
toc

XNplus = XN + zono.cz + sum(abs(zono.Rz),2)*ones(1,N);
XNmoins = XN + zono.cz - sum(abs(zono.Rz),2)*ones(1,N);
cost_rimproved = cost_calcul(Y,U,MPCobj);

%% plot
errors = 0;
outputs = 1;
controls = 1;
observer = 0;
zonoplot = 1;
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
p = [cost_matlab;cost_bob];
bar(p,'stacked')
legend("Pursuit error","Command power","Delta Command")
title("Cost of the trajectory")
xlabel("Matlab      vs      Bob")
ylabel("Total cost")
