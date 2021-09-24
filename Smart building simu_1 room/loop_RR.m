% Exemple of Sun paper 2020 "Resilient MPC"
%% Parameters 
addpath(genpath('..'))
%systeme
init_scenar_ext_RR;
Init_controller_RR;
%% loop 
T = 0:Ts/3600:144;
DoS = ones(length(T),1);
DoS(11:18)= zeros(8,1);
DoS(21:28)= zeros(8,1);
DoS(31:37)= zeros(7,1);
DoS(39) = 0;
DoS(1:2)= zeros(2,1);
x0 = int.x0;
xn_1 = x0;

Ref = [0;0];
Tmax = length(T);
X = zeros(7,Tmax);
Xn = zeros(7,Tmax);
X(:,1) = x0;
Xn(:,1) = x0;
Xn2= zeros(7,Tmax);
U = zeros(1,Tmax);
U_mpc = zeros(1,Tmax);
flag = 0;
u_prec = 0;
mpc_data = struct;
mpc_data.x0 = int.x0;
mpc_data.u0 = int.u0;
mpc_data.Vm = int.vm;
MPCobj.zono.cu_mpc = zono.cu_mpc2;
MPCobj.zono.Ru_mpc = zono.Ru_mpc2;

for i = 1:Tmax-1
%     if in(ZA_opt,X(:,i)-Xn(:,i))
%         MPCobj.zono.cu_mpc = zono.cu_mpc;
%         MPCobj.zono.Ru_mpc = zono.Ru_mpc;
%     else
%         MPCobj.zono.cu_mpc = zono.cu_mpc2;
%         MPCobj.zono.Ru_mpc = zono.Ru_mpc2;
%     end
    if i==1
        [u,mpc_data,xn] = MPC_Robuste(MPCobj,X(:,i),Ref,Vm.Data(:,:,i:end)-Vmd2(:,1),U(:,1),mpc_data,1);
    else
        [u,mpc_data,xn] = MPC_Robuste(MPCobj,X(:,i),Ref,Vm.Data(:,:,i:end)-Vmd2(:,1),U(:,i-1),mpc_data,1);
    end
    if DoS(i)
        U(i) = u(1);
        Ukj = u;
        flag =1;
        nbDoS = 0;
        mpc_data.u_prec = mpc_data.u_prec;
    else
        if flag==0
            U(i)=0;
            mpc_data.u_prec = 0;
        else
            nbDoS = nbDoS+1;
            U(i)=Ukj(nbDoS+1);
            mpc_data.u_prec = Ukj(nbDoS+1);
        end
    end
    
    X(:,i+1) = A*X(:,i) + B*U(i) + D*VD.Data(:,:,i+1);
    Xn2(:,i) = mpc_data.Xnominal;
    Xn(:,i+1) = A*mpc_data.Xnominal + B*mpc_data.u_prec + D(:,2:3)*Vm.Data(:,:,i+1);
    U_mpc(:,i) = mpc_data.u_prec;
    disp(i)
end

Xn(:,N) = Xn(:,N-1);
XN = Xn;
cz = ZA_opt.Z(:,1);
Rz = ZA_opt.Z(:,2:end);
XNplus = Xn + cz + abs(Rz)*ones(size(Rz,2),Tmax);
XNmoins = Xn + cz - abs(Rz)*ones(size(Rz,2),Tmax);

%% plot
dU = U;
dU(end) = dU(end-1);
U_mpc(end) = U_mpc(end-1);
Y = X;
Ref = mpc_data.x0*ones(1,Tmax);
errors = 0;
outputs = 1;
controls = 1;
observer = 0;
zonoplot = 0;
scenar_plot = 0;
plot_simu_script;

%%
figure(2)
plot(U)
grid on
hold on 
plot(U_mpc-U)
%%

figure()
bar(T,1-DoS(1:Tmax))
xlabel("Time(h)")
ylabel("\nu_k")
title("DoS sequence")
grid on
