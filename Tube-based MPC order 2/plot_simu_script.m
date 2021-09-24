
%% plot erreur de poursuite :
if errors
    figure('Name',"error")
    grid on
    for zone = 1:2
        subplot(1,2,zone)
        grid on
        plot(T/3600,Y(zone,:)-ref*ones(1,N),'b')
        title(['Ecart de température réelle zone ',num2str(zone)])
        xlabel('Temps (h)')
        ylabel('Température (°c)')
    end
end
%% plots 1 Y
if outputs
    figure('Name',"Outputs")
    grid on
    for zone = 1:2
        subplot(1,2,zone)
        plot(T/3600,Y(zone,:),'b')
        hold on
        grid on
        %plot(T/3600,dY(zone,:),'b')
        plot(T/3600,XN(zone,:),'r')
        plot(T/3600,XNplus(zone,:),'m')
        plot(T/3600,XNmoins(zone,:),'m')
        plot(T/3600,ref(zone)*ones(1,N),'g')
        legend("modèle réel","modele nominal","borne supérieure","borne inférieure","référence")
        title(['Output ',num2str(zone)])
    end
end

%% plots 2 Command
if controls
    figure('Name',"Control")
    grid on
    for input = 1
        %subplot(2,3,input)
        plot(T/3600,U(input,:),'b')
        hold on
        grid on
        plot(T/3600,dU(input,:),'r')
        legend("Commande envoyé","Commande appliquée")
        title(['Commande',num2str(input)])
        xlabel('Temps (h)')
        ylabel('W/m²')
    end
end
% %% Plot 3 condition extérieurs
% figure('Name',"Disturbance")
% subplot(3,1,1)
% plot(T/3600,V(4,:),T/3600,Vmd(1,:))
% legend("réelle","prédite")
% title('Tempéraure extérieure')
% 
% ylabel('température(°)')
% 
% subplot(3,1,2)
% plot(T/3600,V_RS_E,T/3600,Vmd(2,:))
% hold on
% plot(T/3600,V_RS_W,T/3600,Vmd(3,:))
% legend("Facade est réel","Facade est prédit","Facade ouest réel","Facade ouest prédit")
% title('Ensoleillement facade (W/m²)')
% 
% 
% subplot(3,1,3)
% plot(T/3600,V_IG_bath,T/3600,V_IG_bedroom,T/3600,V_IG_kitchen)
% legend("Salle de bain","Chambres","Cuisines")
% title('Production de chaleur interne (W/m²)')
% xlabel('Temps (h)')
%% observer
if observer
    figure('Name',"Observer")
    grid on
    C = Plant.C;
    Yo = C*Xo;
    for zone = 1:2
        subplot(4,2,zone)
        grid on
        plot(T/3600,Y(zone,:),'b')
        hold on
        plot(T/3600,Yo(zone,:),'r')
        plot(T/3600,ref*ones(1,N),'g')
        legend("Température réelle","température observée","référence")
        title(['Température zone ',num2str(zone)])
        xlabel('Temps (h)')
        ylabel('Température (°c)')
    end
end

%% zono

if zonoplot
   figure('Name',"Zono")
   grid on
   for i =1:N
    drawZ2D(zono.cz+XN(:,i),zono.Rz,'line,b');
    hold on
   end
   legend("Z set")
   plot(XN(1,:),XN(2,:),'b*-',"DisplayName","Nominal")
   plot(Y(1,:),Y(2,:),'ro-',"DisplayName","Real")
   legend()
end