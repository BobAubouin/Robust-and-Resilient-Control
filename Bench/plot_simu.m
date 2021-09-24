%% Condition météo extérieur
t = 0:0.25:72; %temps en heure
V_Text = 14+7.5*cos(2*pi*(t-17)/24);%température extérieure
RS_tot = max(0,900*cos(2*pi*(t-14.5)/24)); %Radiations solaire maximale
id_E = (mod(t,24)<=12)&(mod(t,24)>=8); %index ou les soleil est à l'est (8-12h)
V_RS_E = double(id_E);
V_RS_E(id_E) = 0.8*RS_tot(id_E); %ensoleilllment facade est avec un petit coeff parce qu'on ne tape pas droit
id_W = (mod(t,24)>=16)&(mod(t,24)<=19); %index ou les soleil est à l'ouest (16-19h)
V_RS_W = double(id_W);
V_RS_W(id_W) = 0.8*RS_tot(id_W); %ensoleilllment facade ouest avec un petit coeff parce qu'on ne tape pas droit

%% Condition intérieur

id_bath = (mod(t,24)>=19)&(mod(t,24)<=20);
V_IG_bath = 5*id_bath;

id_bedroom = (mod(t,24)>=21)&(mod(t,24)<=9);
V_IG_bedroom = 5*id_bedroom;

id_kitchen = (mod(t,24)>=19)&(mod(t,24)<=20);
V_IG_kitchen = 3*id_kitchen;

% V = zeros(6,length(t));
% V(4,:) = V_Text;
% V(5,:) = V_RS_E;
V = [V_IG_bath;
     V_IG_bedroom;
     V_IG_kitchen;
     V_Text;
     V_RS_E;
     V_RS_W];
 
Vm = timeseries(V(4:6,:),t*3600);
VD = timeseries(V,t*3600);

%% sim
x0 = 18*ones(size(A_k,1),1);
rx0 = T(1:ordre,:)*x0;
sim('simu_mpc_coloc')


%% plots 1
for zone = 1:8
    subplot(4,2,zone)
    plot(out.tout/3600,out.Y(:,zone),out.tout/3600,20*ones(length(out.tout),1))
    legend("Température réelle","référence")
    title(['Température zone ',num2str(zone)])
    xlabel('Temps (h)')
    ylabel('Température (°c)')
end

figure()
for input = 1:5
    subplot(2,3,input)
    plot(out.tout/3600,out.U(:,input))
    title(['Commande ',num2str(input)])
    xlabel('Temps (h)')
    ylabel('W/m²')
end

%% plots 2
Ylin = reshape(out.Y,[8,289]);
for zone = 1:8
    subplot(4,2,zone)
    plot(out.tout/3600,Ylin(zone,:),out.tout/3600,20*ones(length(out.tout),1))
    legend("Température réelle","référence")
    title(['Température zone ',num2str(zone)])
    xlabel('Temps (h)')
    ylabel('Température (°c)')
end

figure()
for input = 1:5
    subplot(2,3,input)
    plot(out.tout/3600,out.U(:,input))
    title(['Commande ',num2str(zone)])
    xlabel('Temps (h)')
    ylabel('W/m²')
end
%% compare

for zone = 1:8
    subplot(4,2,zone)
    plot(out.tout/3600,out.Y(:,zone),out.tout/3600,20*ones(length(out.tout),1),out.tout/3600,Ylin(zone,:))
    legend("Température réelle","référence","température sur modèle linéaireé")
    title(['Température zone ',num2str(zone)])
    xlabel('Temps (h)')
    ylabel('Température (°c)')
end