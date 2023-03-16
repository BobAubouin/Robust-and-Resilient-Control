% Exemple of Sun paper 2020 "Resilient MPC"
%% Parameters 
addpath(genpath('..'))
%systeme
init_loop_2

%% loop 
Tmax = 30;
Tattack = [1,2,3,4,11,12,13,14,21,23,24,26];

DoS = ones(1,Tmax);
rng(2);
W = (2*rand(1,Tmax)-1)/10;
for i = 1:length(Tattack)
    DoS(Tattack(i))=0;
end
for i = 1:Tmax
    if  i<10
        W(i)=-1/10;
    elseif i<25
        W(i)=1/10;
    end
end


x0 = [1.4;-3];
%xn_1 = x0;
Ref = [0;0];

X = zeros(2,Tmax);
Xn = zeros(2,Tmax);
X(:,1) = x0;
Xn(:,1) = x0;
U = zeros(1,Tmax);
U_mpc = zeros(1,Tmax);
flag = 0;
u_prec = 0;
mpc_data = struct;
mpc_data.u0 = 0;
mpc_data.x0 = [0;0];
mpc_data.Vm = 0;

for i = 1:Tmax-1
    if i==1
        [u,mpc_data,xn] = MPC_Robuste_poly(MPCobj,X(:,i),Ref,[],U(:,1),mpc_data,1);
    else
        [u,mpc_data,xn] = MPC_Robuste_poly(MPCobj,X(:,i),Ref,[],U(:,i-1),mpc_data,1);
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
    
    X(:,i+1) = A*X(:,i) + B*U(i) + D*W(i);
    Xn(:,i+1) = A*mpc_data.Xnominal + B*mpc_data.u_prec;
    U_mpc(:,i) = mpc_data.u_prec;
    disp(i)
end

figure()

plot(Xn_4,dim,'g');
hold on
plot(X(1,:),X(2,:),'bx-')
grid on
plot(Xn(1,:),Xn(2,:),'ro--')
%plot(ZX,[1 2]);
for i=1:Tmax
    plot(Xn(:,i)+ZA_opt,[1 2],'r');
end

legend("Initial feasible set (4 DoS)",'Real system',"Nominal system","$N$-RPI tube",'Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
xlabel('$x_1$','Interpreter','latex')



figure()

plot(X(1,:)-Xn(1,:),X(2,:)-Xn(2,:),'bx-')
grid on
hold on
plot(ZA_opt,[1 2],'r');
plot(ZB_opt,[1 2],'r--');
plot(RPI,[1 2],'g');
ylabel('$z_2$','Interpreter','latex')
xlabel('$z_1$','Interpreter','latex')
legend("State trajectories ($z_k$)","N-RPI set ($\bf{Z}$)","Bound of the trajectory ($\bf{Z^+}$)","mRPI set ($\bf{P}$)",'Interpreter','latex')
% axes('position',[.65 .175 .25 .25])
% box on % put box around new pair of axes
% plot(X(1,:)-Xn(1,:),X(2,:)-Xn(2,:),'bx-')
% grid on
% hold on
% plot(ZA_opt,[1 2],'r');
% plot(RPI,[1 2],'g');
%axis tight

%%
figure()

plot(U,'b')
grid on
hold on 
plot(U-U_mpc,'r')
legend("Total control ($u_k$)","Proportionnal input($K(x_k-\bar{x}_k)$)",'Interpreter','latex')
xlabel("$k$",'Interpreter','latex')
%%
figure()
plot(0:Tmax-1,W(1:Tmax))
xlabel("$k$",'Interpreter','latex')
ylabel("$w_k$",'Interpreter','latex')
grid on
%%

figure()
bar(0:Tmax-1,1-DoS(1:Tmax))
hold on
%bar(0:Tmax-1,DoS(1:Tmax))

xlabel("$k$",'Interpreter','latex')
ylabel("$1-\nu_k$",'Interpreter','latex')
%title("DoS sequence")

%% states

figure()

Xrien = X;
Xpert = X;
plot(0:Tmax-1,Xrien(1,:),'b:')
hold on
plot(0:Tmax-1,Xrien(2,:),'r:')

plot(0:Tmax-1,Xpert(1,:),'b--')
plot(0:Tmax-1,Xpert(2,:),'r--')

plot(0:Tmax-1,X(1,:),'b')
plot(0:Tmax-1,X(2,:),'r')
xlabel("$k$",'Interpreter','latex')
grid on
legend("$x_1$ Without DoS, $w_k=0$","$x_2$ Without DoS, $w_k=0$","$x_1$ Without DoS, $w_k \ne 0$ ", "$x_2$ Without DoS, $w_k \ne 0$", "$x_1$ With DoS, $w_k \ne 0$ ", "$x_2$ With DoS, $w_k \ne 0$",'Interpreter','latex')