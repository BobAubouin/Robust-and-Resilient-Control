%% Simulation test

%% Condition météo extérieur
addpath(genpath('..'))
init_scenar_ext
 

%% get full model

load("discrete_model_coloc") % data from BRCM

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

% Xmean = mean(X,2);
% Vmean = mean(V,2);
load('datalin.mat')

Blin = bilin2lin(Bu,Bv,Bxu,Bvu,Xmean,Vmean);
Blin_k = bilin2lin(Bu_k,Bv_k,Bxu_k,Bvu_k,Xmean,Vmean);
Blin_k(:,1) = Blin_k(:,1)/4;
Blin_k(:,2) = Blin_k(:,2)/4;
syslin = ss(A_k,Blin_k,C_k,[],60*15,'StateName',B.building_model.identifiers.x);


ordre = 12;

[sysb,g,Tt,Ti] = balreal(syslin);


drsys_tc = modred(sysb,ordre+1:size(A,1),'Truncate');


drsys_tc_12 = drsys_tc;

save("reduced_model_coloc",'drsys_tc_12','Tt')
save("full_li_model_coloc",'syslin')
%% Simulate 

load("reduced_model_coloc",'drsys_tc_12','Tt')
load("full_li_model_coloc",'syslin')
drsys_tc = drsys_tc_12;
Ts = 15*60; %discretization toute les 15 minutes
t = 0:Ts/3600:72;

x0 = 22*ones(size(A_k,1),1);
[Y_full,t2] = simbilin(A_k,Bu_k,Bv_k,Bxu_k,Bvu_k,C_k,t,U,V,x0);

u = [U;V];
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
% figure()
% plot(t,Y_full(zone,1:end-1)-y_lin(:,zone)')
