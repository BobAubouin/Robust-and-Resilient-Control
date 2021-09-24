% Exemple of Sun paper 2020 "Resilient MPC"
%% Parameters 
addpath(genpath('..'))
%systeme
init_loop

%% loop 
DoS = ones(1,30);
DoS(1) = 0;
DoS(3) = 0;
%DoS(5) = 0;
%DoS(7:10) = 0;
% DoS(8) = 0;
% DoS(9) = 1;
% DoS(10) = 1;
 DoS(12) = 0;
 DoS(13) = 0;
% DoS(18) = 0;
% DoS(19) = 0;
% DoS(20) = 0;
% DoS(21) = 1;
% 
% DoS(23) = 0;
% DoS(27) = 0;


x0 = [7;7];
xn_1 = x0;
Ref = [0;0];
Tmax = 30;
X = zeros(2,Tmax);
Xn = zeros(2,Tmax);
X(:,1) = x0;
Xn(:,1) = x0;
U = zeros(1,Tmax);
U_mpc = zeros(1,Tmax);
flag = 0;
u_prec = 0;
mpc_data = struct;
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
        [u,mpc_data,xn] = MPC_Robuste(MPCobj,X(:,i),Ref,[],U(:,1),mpc_data,1);
    else
        [u,mpc_data,xn] = MPC_Robuste(MPCobj,X(:,i),Ref,[],U(:,i-1),mpc_data,1);
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
    
    X(:,i+1) = A*X(:,i) + B*U(i) + D*(Zd.center-Zd.generators*(2*rand(2,1)-1));
    Xn(:,i+1) = A*mpc_data.Xnominal + B*mpc_data.u_prec;
    U_mpc(:,i) = mpc_data.u_prec;
    disp(i)
end

figure()

plot(Xn_2,dim,'FaceColor',[0 .6 0.8],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
title("State trajectories with Z set")
hold on
plot(X(1,:),X(2,:),'bo-')
grid on
plot(Xn(1,:),Xn(2,:),'ro--')
for i=1:Tmax
    plot(Xn(:,i)+ZA_opt,[1 2],'r--');
end
%'$ \bf{X}_{10} $','$ \bf{X}_{10}^1 $','$ \bf{X}_{10}^2 $',
legend("Initial feasible set (2 DoS)",'Real system',"Nominal system and $\mu$RPI",'Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
xlabel('$x_1$','Interpreter','latex')

%%
figure()

plot(U)
grid on
hold on 
plot(U_mpc-U)
legend("Total control","Proportionnal input")
title("Control input")
%%
figure()

plot(X(1,:))
hold on
grid on
plot(X(2,:))
legend('x_1','x_2')
title("States trajectories")

figure()
bar(1-DoS(1:Tmax))
xlabel("k")
ylabel("\nu_k")
title("DoS sequence")
