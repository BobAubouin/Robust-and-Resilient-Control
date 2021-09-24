
init_scenar_ext;
init_plant;
%Init_full_controller;
Init_controller;
init_attacker_policy;
%% init
tic
choix_scenario =6; %1 = pas d'attaque, 2 = DoS, 3 = biais, 4 = upper sat, 5 = replay, 6 false prediction

ordre = size(MPCobj.Model.Plant.A,2);
x0 = [   19.6733
       16.6278
       14.5857
       14.0955
       14.0955
       14.0955
       16.0827];
ref = 20*ones(1,1);
T = 0:0.25:72;
N = length(T);
U = zeros(size_input,N);
dU = zeros(size_input,N);
Y = zeros(size_outputs,N);
dY = zeros(size_outputs,N);
Y(:,1) = x0(1);
dY(:,1) = x0(1);
X = zeros(size_state,N);
X(:,1) = x0;
Xo = zeros(ordre,N);
Xo(:,1) = x0;
xmpc = mpcstate(MPCobj,Xo(:,1));

BIY = variable0.BIY;
BIU = variable0.BIU;
Replay_data = variable0.RAYsave;


for i = 1:N-1
   t = T(i);
   %[X_k1,y] = Bilin(X(:,i),dU(:,i),VD.Data(:,:,i),A_k,Bu_k,Bv_k,Bxu_k,Bvu_k,C_k);
   X_k1 = A_k*X(:,i) + [Bu_k,Bv_k]*[dU(:,i);VD.Data(:,:,i)];
   y = C_k*X_k1;
   X(:,i+1) = X_k1;
   Y(:,i+1) = y;
   
   [dy,BIY,Replay_data,pred] = attacker_Y(y,Vm.Data(:,:,i+1:end),t,BIY,Replay_data,choix_scenario,Data,Ts);
   dY(:,i+1) = dy;
   
   xmpc.Plant = X_k1;
   [u] = control_adaptive(MPCobj,xmpc,dy,ref,pred);
   U(:,i+1) = u;
   Xo(:,i+1) = xmpc.Plant;
    
   [du,BIU] = attacker_U(u,dU(:,i),t,BIU,choix_scenario,Data,Ts);
    dU(:,i+1) = du;
end
toc
%% plot

Y(2) = (Y(1)+Y(3))/2;
plot(T,Y);
hold on
%plot(T,X(5,:));
xlabel("Time (h)")
ylabel("Interior Temperature (°)")
%legend("Interior temperature","Wall temperature")

figure()
plot(T(2:end),U(2,2:end)*25);
hold on
%plot(T(2:end),dU(2,2:end)*25);
xlabel("Time (h)")
ylabel("Heat power (W)")
%legend("Calculated input","Applied input")

%%

errors = 0;
outputs = 1;
controls = 0;
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
