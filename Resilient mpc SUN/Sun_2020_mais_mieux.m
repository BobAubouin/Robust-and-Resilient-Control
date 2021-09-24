% Correction of Sun paper 2020 "Resilient MPC"

% This script correct the numerical example present in the article of SUN 2020.


%% Parameters 
addpath(genpath('..'))
%systeme
A = [1   -1.2; 1.2 1.1];
B = [1 ; 0.5]; %NOK

% MPC controller
Q = eye(2);
R = 0.1;
P = [1.9385, 1.7088 ; 1.7088, 5.0552];
K = [-1.5718,- 0.2611];
N = 10;
M = 2;

[h1,h2,lambda1,lambda2,ok] = get_param_h_l(A,B,K,P,N,M,1);



% set 

Zu = zonotope([0 1]);
Zx = zonotope([0 10 0;
               -3 0 5]);
Zw = zonotope([0 0.1;
               0 0.1]);
           
           
alpha = TerminalSet(K,P,Zu,Zx);
hlambda_m = max(h1*lambda1,h2*lambda2);
coeff = zeros(N,1);
for i=1:N
    if i <= M
        coeff(i) = (hlambda_m)^i;
    else 
        coeff(i) = (h1^M*h2^min((M+1),i-M)*lambda1^M*lambda2^(i-M));
    end
end
h_star = max(coeff);
alpha_bar = alpha/h_star^2;

Eterm = ellipsoid(inv(P)*alpha_bar);
Eterm_init = ellipsoid(inv(P)*alpha);
Eterm_sun = ellipsoid(inv(P)*0.1);
Zterm = zonotope(Eterm,5,'i:norm');
Zterm_sun = zonotope(Eterm_sun,5,'i:norm');
           
Xn = InverseReach_1(Zterm,A,B,Zu,N,0);
Xn_1 = InverseReach_1(Zterm,A,B,Zu,N,1);
Xn_2 = InverseReach_1(Zterm,A,B,Zu,N,2);

Xn_sun = InverseReach_1(Zterm_sun,A,B,Zu,N,0);
Xn_1_sun = InverseReach_1(Zterm_sun,A,B,Zu,N,1);
Xn_2_sun = InverseReach_1(Zterm_sun,A,B,Zu,N,2);

init_controller_sun_2020
%% plot set
figure(1)
%plot(Eterm_init,[1 2],'FaceColor',[0 0.6 0],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
hold on
axis equal
grid on
plot(Eterm,[1 2],'FaceColor',[.6 0 0],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
% xlim([-1.5,1.5]);
% ylim([-1.5,1.5]);
%plot(Xn_sun,[1 2],'FaceColor',[0.6 0 0],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
%plot(Xn,[1 2],'FaceColor',[0 .6 0],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
% plot(Xn_1_sun,[1 2],'FaceColor',[0.6 .6 0.6],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
%plot(Xn_1,[1 2],'FaceColor',[0 .6 0.6],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
% plot(Xn_2_sun,[1 2],'FaceColor',[0.6 0 .6],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
%plot(Xn_2,[1 2],'FaceColor',[0 0 .6],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
%plot(Zterm_sun,[1 2],'FaceColor',[ 0 .6 0],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
%plot(Zterm,[1 2],'FaceColor',[.6 0 0],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
%% loop 
DoS = ones(1,30);
DoS(1) = 0;
DoS(3) = 0;
%DoS(5) = 0;
%DoS(7:10) = 0;
% DoS(8) = 0;
% DoS(9) = 1;
% DoS(10) = 1;
 DoS(11) = 1;
 DoS(12) = 1;
% DoS(18) = 0;
% DoS(19) = 0;
% DoS(20) = 0;
% DoS(21) = 1;
% 
% DoS(23) = 0;
% DoS(27) = 0;

%x0 = [-0.25;0.05];
x0 = [-0.046;0.014];
%x0 = [0.044;-0.022];
%x0 = [0.14;-0.08];
%x0 = [-1,0.5];
Tmax = 10;
X = zeros(2,Tmax);
X(:,1) = x0;
U = zeros(1,Tmax);
Unsat = zeros(1,Tmax);
flag = 0;
u_prec = 0;
mpc_data = struct;
for i = 1:Tmax-1
    if in(Eterm,X(:,i)) || flag
        flag = 1;
        u = max(min(K*X(:,i)*DoS(i),1),-1);
        unsat = K*X(:,i)*DoS(i);
        disp(i)
    elseif i>2
        if i>1
            u_prec = U(i-1);
        end
        [u,mpc_data] = MPC_zono(MPCobj,X(:,i),[],u_prec,mpc_data,zono); 
        unsat = u;
        %disp(u)
    else
        u = 0;
        unsat = 0;
    end
    %u = K*X(:,i)*DoS(i);
    U(i) = u;
    Unsat(i) = unsat;
    X(:,i+1) = A*X(:,i) + B*u;% + [1;1]*(rand()*2-1)*0.1;
    
end

plot(X(1,:),X(2,:),'bo-')
%'$ \bf{X}_{10} $','$ \bf{X}_{10}^1 $','$ \bf{X}_{10}^2 $',
legend('$ \bf{X}_{f} $',"State trajectory",'Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
xlabel('$x_1$','Interpreter','latex')

figure(2)
plot(U)
grid on
hold on 
plot(Unsat)
legend("Saturated control input","Calculated control input")
figure(3)
plot(X(1,:))
hold on
grid on
plot(X(2,:))
legend('x_1','x_2')

figure()
bar(1-DoS(1:Tmax))
xlabel("k")
ylabel("\nu_k")
title("DoS sequence")
