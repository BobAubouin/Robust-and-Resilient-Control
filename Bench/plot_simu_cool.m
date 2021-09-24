%% plot erreur de poursuite :
ref = out.tout*0+20;
N = length(out.tout);
figure('Name',"error")
for zone = 1:8
    subplot(4,2,zone)
    plot(out.tout/3600,out.Y(:,zone)-ref,'b')
    title(['Ecart de température réelle zone ',num2str(zone)])
    xlabel('Temps (h)')
    ylabel('Température (°c)')
end

%% plots 1 Y
figure('Name',"Outputs")
for zone = 1:8
    subplot(4,2,zone)
    plot(out.tout/3600,out.Y(:,zone),'b')
    hold on
    plot(out.tout/3600,out.dY(:,zone),'r')
    plot(out.tout/3600,ref,'g')
    legend("Température réelle","température reçu","référence",'Location','southeast')
    title(['Température zone ',num2str(zone)])
    xlabel('Temps (h)')
    ylabel('Température (°c)')
end
%% plots 2 Command
figure('Name',"Control")
for input = 1:5
    subplot(2,3,input)
    plot(out.tout/3600,out.U(:,input),'b')
    hold on
    plot(out.tout/3600,out.dU(:,input),'r')
    legend("Commande envoyé","Commande appliquée")
    title(['Commande',num2str(input)])
    xlabel('Temps (h)')
    ylabel('W/m²')
end

%% Plot 3 condition extérieurs
figure('Name',"Disturbance")
subplot(3,1,1)
plot(out.tout/3600,V(4,:),out.tout/3600,Vmd(1,:))
legend("réelle","prédite")
title('Tempéraure extérieure')

ylabel('température(°)')

subplot(3,1,2)
plot(out.tout/3600,V_RS_E,out.tout/3600,Vmd(2,:))
hold on
plot(out.tout/3600,V_RS_W,out.tout/3600,Vmd(3,:))
legend("Facade est réel","Facade est prédit","Facade ouest réel","Facade ouest prédit")
title('Ensoleillement facade (W/m²)')


subplot(3,1,3)
plot(out.tout/3600,V_IG_bath,out.tout/3600,V_IG_bedroom,out.tout/3600,V_IG_kitchen)
legend("Salle de bain","Chambres","Cuisines")
title('Production de chaleur interne (W/m²)')
xlabel('Temps (h)')
