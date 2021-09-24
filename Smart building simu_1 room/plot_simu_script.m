
%% plot erreur de poursuite :
if errors
    figure('Name',"error")
    grid on
    for zone = 1:2
        subplot(1,2,zone)
        grid on
        plot(T,Y(zone,:)-Ref*ones(1,N),'b')
        title(['Ecart de température réelle zone ',num2str(zone)])
        xlabel('Temps (h)')
        ylabel('Température (°c)')
    end
end
%% plots 1 Y
if outputs
    figure('Name',"Outputs")
    grid on
    for zone = 1:7
        subplot(2,4,zone)
        
        plot(T,Ref(zone,:),'g')
        hold on
        grid on
        %plot(T,dY(zone,:),'b')
        
        plot(T,XNplus(zone,:),'m')
        plot(T,XNmoins(zone,:),'m')
        plot(T,XN(zone,:),'r-')
        plot(T,Y(zone,:),'b')
        %title(B.building_model.identifiers.x(zone))
    end
    legend("modele nominal","borne supérieure","borne inférieure","référence","modèle réel")
end

%% ZOne temp
figure()
plot(T,Ref(1,:),'g')
hold on
grid on
plot(T,XNplus(1,:),'m')
plot(T,XNmoins(1,:),'m')
plot(T,XN(1,:),'r-')
plot(T,Y(1,:),'b')
xlabel("Time (h)")
ylabel("Temperature (°)")
legend("Target","Upper bound","Lower bound","Nominal model","Real model")

%% plots 2 Command
if controls
    figure('Name',"Control")
    grid on
    for input = 1:1
        subplot(1,1,input)
        plot(T,dU(input,:)*25,'b')
        hold on
        grid on
        plot(T,U_mpc(input,:)*25,'r--')
        legend("Total input","MPC part of the input")
        %title(['Heat power',num2str(input)])
        xlabel('Times (h)')
        ylabel('Heat power (W)')
    end
end

%% Plot 3 condition extérieurs
if scenar_plot
    figure('Name',"Disturbance")
    subplot(3,1,1)
    plot(T,V(2,:),T,Vmd(1,:))
    legend("Real","Predicted")
    title('Exterior Temperature(°)')

    subplot(3,1,2)
    plot(T,V(3,:),T,Vmd(2,:))
    legend("Real","Predicted")
    title('Solar Radiation (W/m²)')


    subplot(3,1,3)
    plot(T,V_IG)
    title('Internal Heat Gain (W/m²)')
    xlabel('Time (h)')
end

%% observer
if observer
    figure('Name',"Observer")
    grid on
    C = Plant.C;
    Yo = C*Xo;
    for zone = 1:2
        subplot(4,2,zone)
        grid on
        plot(T,Y(zone,:),'b')
        hold on
        plot(T,Yo(zone,:),'r')
        plot(T,Ref*ones(1,N),'g')
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