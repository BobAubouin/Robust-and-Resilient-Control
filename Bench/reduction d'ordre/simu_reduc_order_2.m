%% Simulation test

%% Condition météo extérieur
addpath(genpath('..'))
init_scenar_ext
 

%% get full model

load("discrete_model_coloc_3") % data from BRCM

%full discrete systeme
A_k = B.building_model.discrete_time_model.A; % same for Bu,Bv,Bvu,Bxu
Bu_k = B.building_model.discrete_time_model.Bu;
Bv_k = B.building_model.discrete_time_model.Bv;
Bvu_k = B.building_model.discrete_time_model.Bvu;
Bxu_k = B.building_model.discrete_time_model.Bxu;
C_k = diag((1:size(A_k,1))<=8);
C_k = double(C_k(1:8,:));
%% reduced model

A = B.building_model.continuous_time_model.A; % same for Bu,Bv,Bvu,Bxu
Bu = B.building_model.continuous_time_model.Bu;
Bv = B.building_model.continuous_time_model.Bv;
Bvu = B.building_model.continuous_time_model.Bvu;
Bxu = B.building_model.continuous_time_model.Bxu;

%Xmean = mean(X,2);
% Vmean = mean(V,2);
% Umean = U(:,end);
% Ucst = repmat(Umean,1,size(V,2));
% Ucst = zeros(5,289);

load('datalin.mat')
Umean = [Umean;5;5;5];

Blin = bilin2lin(Bu,Bv,Bxu,Bvu,Xmean,Vmean);
Blin_k = bilin2lin(Bu_k,Bv_k,Bxu_k,Bvu_k,Xmean,Vmean);
Blin_k(:,1) = Blin_k(:,1)/4;
Blin_k(:,2) = Blin_k(:,2)/4;
% Blin_k(:,10) = Blin_k(:,10)*4;
% Blin_k(:,11) = Blin_k(:,11)*4;
syslin = ss(A_k,[Bu_k,Bv_k],C_k,[],60*15,'StateName',B.building_model.identifiers.x);


ordre = 12;

[sysb,g,Tt,Ti] = balreal(syslin);


drsys_tc = modred(sysb,ordre+1:size(A,1),'Truncate');


drsys_tc_12 = drsys_tc;

save("reduced_model_coloc",'drsys_tc_12','Tt')
save("full_li_model_coloc",'syslin')
%% Simulate with command echelon

drsys_tc = drsys_tc_12;
Ts = 15*60; %discretization toute les 15 minutes
t = 0:Ts/3600:96;
if exist('X','var') == 1
    x0 = X(:,end);
else
    x0 = 20*ones(133,1);
end
Uech = repmat([0;0;4;4;4;40;4;4;4],1,length(t));
Uech(6,100:end) = 0*ones(1,length(t)-99);
V2=[0;0;0;15;0;0];
V2 = repmat(V2(:,1),1,length(t));

%x0 = X(:,end);
[Y_full,t2,X] = simbilin(A_k,Bu_k,Bv_k,Bxu_k,Bvu_k,C_k,t,Uech,V2,x0);

u = [Uech;V2];
if exist('X2','var') == 1
    x0 = X2(end,:);
    rx0 = Tt(1:ordre,:)*x0';
else
    x0 = 20*ones(133,1);
    rx0 = Tt(1:ordre,:)*x0;
end
[y_lin,~,X2] = lsim(syslin,u,t*3600,x0);



[y_reduce_tc,~,X3] = lsim(drsys_tc,u,t*60*60,rx0);


for zone = 1:8
    subplot(4,2,zone)
    plot(t,y_lin(:,zone),t,y_reduce_tc(:,zone))
    legend("modele linéaire","modele reduit")
    title(['Température zone ',num2str(zone)])
    xlabel('Temps (h)')
    ylabel('Température (°c)')
end
figure()
zone = 7;
plot(t,y_lin(:,zone),t,y_reduce_tc(:,zone))
legend("linear model","reduced model")
%title(['Température zone ',num2str(zone)])
xlabel('Time (h)')
ylabel('Temperature (°c)')
% plot(t,Y_full(zone,1:end-1)-y_lin(:,zone)')

%% Simulate with distubance echelon

drsys_tc = drsys_tc_12;
Ts = 15*60; %discretization toute les 15 minutes
t = 0:Ts/3600:5000;

Uech = repmat([1;0.5;3.2;3.2;3.2;3.2;3.2;3.2;3.2],1,length(t));
%Uech(3,10001:11000) = 50*ones(1,1000);
V2 = repmat(V2(:,1),1,length(t));
V2(5,10001:11000) = V2(5,10001:11000)-1;

x0 = Xmean;
%x0 = X(:,end);
[Y_full,t2,X] = simbilin(A_k,Bu_k,Bv_k,Bxu_k,Bvu_k,C_k,t,Uech,V2,x0);

u = [Uech;V2];
y_lin = lsim(syslin,u,t*3600,x0);

rx0 = Tt(1:ordre,:)*x0;

y_reduce_tc = lsim(drsys_tc,u,t*60*60,rx0);


for zone = 1:8
    subplot(4,2,zone)
    plot(t,Y_full(zone,1:end-1),t,y_lin(:,zone),t,y_reduce_tc(:,zone))
    legend("modele complet","modele linéaire","modele reduit")
    title(['Température zone ',num2str(zone)])
    xlabel('Temps (h)')
    ylabel('Température (°c)')
end
figure()
zone = 7;
plot(t,Y_full(zone,1:end-1),t,y_lin(:,zone),t,y_reduce_tc(:,zone))
legend("modele complet","modele linéaire","modele reduit")
title(['Température zone ',num2str(zone)])
xlabel('Temps (h)')
ylabel('Température (°c)')
% plot(t,Y_full(zone,1:end-1)-y_lin(:,zone)')