%% Condition météo extérieur
Ts = 15*60; %discretization toute les 15 minutes
T = 0:Ts/3600:144; %temps en heure
V_Text = 10+2*cos(2*pi*(T-17)/24);%température extérieure
RS_tot = max(0,900*cos(2*pi*(T-14.5)/24)); %Radiations solaire maximale (W/m²)
id_E = (mod(T,24)<=12)&(mod(T,24)>=8); %index ou les soleil est à l'est (8-12h)
V_RS_E = double(id_E);
V_RS_E(id_E) = 0*RS_tot(id_E); %ensoleilllment facade est avec un petit coeff parce qu'on ne tape pas droit

%% Condition intérieur

id = (mod(T,24)>=18)&(mod(T,24)<=20);
V_IG = 0*id;% + 5*rand(1,length(t));


% V = zeros(6,length(t));
% V(4,:) = V_Text;
% V(5,:) = V_RS_E;
V = [V_IG;
     V_Text+4*(rand(1,length(T))-0.5);
     V_RS_E.*(1+0.4*(rand(1,length(T))*2-1))];
V2 = repmat(mean(V,2),1,size(V,2));



Vmd = [V_Text;
        V_RS_E ];
    
Vm = timeseries(Vmd,T*3600);

Vmd2 = [repmat(mean(V_Text),1,size(V,2));
        repmat(mean(V_RS_E),1,size(V,2))];
    
Vm2 = timeseries(Vmd2,T*3600);


VD = timeseries(V,T*3600);
VD2 = timeseries(V2,T*3600);