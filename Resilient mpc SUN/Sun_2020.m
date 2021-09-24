% Exemple of Sun paper 2020 "Resilient MPC"

% This script reproduce the numerical example present in the article of SUN 2020.
% With x0 = [-0.25,0,05]' we obtain instability which show us a problem in their example.

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

% set 

alpha = 0.1;
Eterm = ellipsoid(inv(P)*alpha);
%E2 = ellipsoid(eye(2)*alpha);
Zterm = zonotope(Eterm,5,'i:norm');
Zu = zonotope([0 1]);
Zx = zonotope([0 10 0;
               -3 0 5]);
Zw = zonotope([0 0.1;
               0 0.1]);
           
Xn = InverseReach_1(Zterm,A,B,Zu,N,0);
Xn_1 = InverseReach_1(Zterm,A,B,Zu,N,1);
Xn_2 = InverseReach_1(Zterm,A,B,Zu,N,2);


init_controller_sun_2020
%% plot set 2
figure(1)
%plot(Eterm,[1 2]);
hold on
axis equal
grid on
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
plot(Xn,[1 2],'FaceColor',[0 .6 0],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
plot(Xn_1,[1 2],'FaceColor',[0 .6 0.6],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
plot(Xn_2,[1 2],'FaceColor',[0 0 .6],'Filled',true,'FaceAlpha',0.2,'EdgeColor','none');
plot(Zterm,[1 2],'FaceColor',[.6 0 0],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none');
%plot(E2,[1 2],'FaceColor',[.6 0.6 0],'Filled',true,'FaceAlpha',0.5,'EdgeColor','none')
%% loop 
DoS = ones(1,30);
DoS(1:1) = 0;
DoS(3) = 0;

x0 = [-0.25;0.05];
X = zeros(2,30);
X(:,1) = x0;
U = zeros(1,30);
flag = 0;
u_prec = 0;
mpc_data = struct;
for i = 1:29
    if in(Zterm,X(:,i)) || flag
        flag = 1;
        u = max(min(K*X(:,i)*DoS(i),1),-1);
        disp(i)
    elseif i>0
        if i>1
            u_prec = U(i-1);
        end
        [u,mpc_data] = MPC_zono(MPCobj,X(:,i),[],u_prec,mpc_data,zono); 
        %disp(u)
    else
        u = 0;
    end
    %u = K*X(:,i)*DoS(i);
    U(i) = u;
    
    X(:,i+1) = A*X(:,i) + B*u;
    
end

plot(X(1,:),X(2,:),'bo-')
legend('$ \bf{X}_{10} $','$ \bf{X}_{10}^1 $','$ \bf{X}_{10}^2 $','$ \bf{X}_{f} $','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
xlabel('$x_1$','Interpreter','latex')

figure(2)
plot(U)
grid on
figure(3)
plot(X(1,:))
hold on
grid on
plot(X(2,:))
legend('$x_1$','$x_2$','Interpreter','latex')



