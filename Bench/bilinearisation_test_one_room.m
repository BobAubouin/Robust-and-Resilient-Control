load('discrete_model_solo.mat')
A_k = B.building_model.discrete_time_model.A; % same for Bu,Bv,Bvu,Bxu
Bu_k = B.building_model.discrete_time_model.Bu;
Bv_k = B.building_model.discrete_time_model.Bv;
Bvu_k = B.building_model.discrete_time_model.Bvu;
Bxu_k = B.building_model.discrete_time_model.Bxu;
C_k = [1 0 0 0 0 0 0];
%%

Vmean = [0; 14; 60];
%
Ts = 15*60;
t = 0:Ts/3600:200;

x0 =[18.65
   16.7116
   23.8350
   15.8039
   15.8039
   15.8039
   16.2742];

U = repmat([0.5;25],1,length(t));
V = repmat(Vmean,1,length(t));
[Y_full,t2,X] = simbilin(A_k,Bu_k,Bv_k,Bxu_k,Bvu_k,C_k,t,U,V,x0);
plot(t2,Y_full)
xlim([0,150])
xlabel("Time (hour)");
ylabel("Temperature of the room (°c)");
grid on
Xmean = X(:,end);

Blin = bilin2lin(Bu_k,Bv_k,Bxu_k,Bvu_k,Xmean,Vmean);

%%
t = 0:Ts/3600:24;

Vmean = [0; 14; 100];
V = repmat(Vmean,1,length(t));



RS_tot = max(0,900*cos(2*pi*(t-14.5)/24)); %Radiations solaire maximale (W/m²)
id_E = (mod(t,24)<=12)&(mod(t,24)>=8); %index ou les soleil est à l'est (8-12h)
V_RS_E = double(id_E);
V_RS_E(id_E) = RS_tot(id_E);
V(3,:) = V_RS_E;
% 
% plot(t,V_RS_E)
% xlabel("Time (hour)");
% ylabel("Solar Radiation (W/m^2)");
% grid on



U_blind = [repmat([1;25],1,length(t)/2+0.5), repmat([1;25],1,length(t)/2-0.5)];
[Y_bili,t2,X] = simbilin(A_k,Bu_k,Bv_k,Bxu_k,Bvu_k,C_k,t,U_blind,V,Xmean);
[Y_li,t2,X] = simlin(A_k,Blin(:,1:2),Bv_k,C_k,t,U_blind,V,Xmean);

plot(t2,Y_bili);
hold on
plot(t2,Y_li,'*');

xlabel("Time (hour)");
ylabel("Temperature of the room (°c)");
grid on
legend("Bilinear model", "Linearized model")

