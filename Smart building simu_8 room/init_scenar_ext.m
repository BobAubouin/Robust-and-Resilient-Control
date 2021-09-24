%% Condition météo extérieur
Ts = 15*60; %discretization toute les 15 minutes
t = 0:Ts/3600:72; %temps en heure
V_Text = -25+5*cos(2*pi*(t-17)/24);%température extérieure
RS_tot = max(0,900*cos(2*pi*(t-14.5)/24)); %Radiations solaire maximale (W/m²)
id_E = (mod(t,24)<=12)&(mod(t,24)>=8); %index ou les soleil est à l'est (8-12h)
V_RS_E = double(id_E);
V_RS_E(id_E) = 0.2*RS_tot(id_E); %ensoleilllment facade est avec un petit coeff parce qu'on ne tape pas droit
id_W = (mod(t,24)>=16)&(mod(t,24)<=19); %index ou les soleil est à l'ouest (16-19h)
V_RS_W = double(id_W);
V_RS_W(id_W) = 0.2*RS_tot(id_W); %ensoleilllment facade ouest avec un petit coeff parce qu'on ne tape pas droit

%% Condition intérieur

id_bath = (mod(t,24)>=18)&(mod(t,24)<=20);
V_IG_bath = 5*id_bath;% + 5*rand(1,length(t));

id_bedroom = (mod(t,24)>=21)|(mod(t,24)<=9);
V_IG_bedroom = 5*id_bedroom;%  + 5*rand(1,length(t));

id_kitchen = (mod(t,24)>=20)&(mod(t,24)<=22);
V_IG_kitchen = 5*id_kitchen;%  + 5*rand(1,length(t));

% V = zeros(6,length(t));
% V(4,:) = V_Text;
% V(5,:) = V_RS_E;
V = [V_IG_bath;
     V_IG_bedroom;
     V_IG_kitchen;
     V_Text+2*(rand(1,length(t))-0.5);
     V_RS_E.*(1+0.2*(rand(1,length(t))*2-1));
     V_RS_W.*(1+0.2*(rand(1,length(t))*2-1))];
 
V2 = repmat(mean(V,2),1,size(V,2));



Vmd = [ repmat(mean(V_IG_bath,2),1,size(V_IG_bath,2));
        repmat(mean(V_IG_bedroom,2),1,size(V_IG_bedroom,2));
        repmat(mean(V_IG_kitchen,2),1,size(V_IG_kitchen,2));
        V_Text;
        V_RS_E;
        V_RS_W];
    
Vm = timeseries(Vmd,t*3600);

Vmd2 = [repmat(mean(V_IG_bath,2),1,size(V_IG_bath,2));
        repmat(mean(V_IG_bedroom,2),1,size(V_IG_bedroom,2));
        repmat(mean(V_IG_kitchen,2),1,size(V_IG_kitchen,2));
        repmat(mean(V_Text),1,size(V,2));
        repmat(mean(V_RS_E),1,size(V,2));
        repmat(mean(V_RS_W),1,size(V,2))];
    
Vm2 = timeseries(Vmd2,t*3600);


VD = timeseries(V,t*3600);
VD2 = timeseries(V2,t*3600);