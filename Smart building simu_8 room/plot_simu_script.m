
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
YN = C*XN; 
Yref = C*Ref;
Name_room = ["Bedroom 1","Bedroom 2","Bedroom 3","Hallway","Bedroom 4","Bathroom 2","Lounge","Bathroom1"];
if outputs
    figure('Name',"Outputs")
    grid on
    for zone = 1:8
        if zone==4
            subplot(4,2,8)
        elseif zone==5
            subplot(4,2,4)
        elseif zone==8
            subplot(4,2,5)
        else
            subplot(4,2,zone)
        end
        plot(T,Yref(zone)*ones(1,N),'g')
        hold on
        grid on
        %plot(T,dY(zone,:),'b')
        plot(T,YN(zone,:),'r')
        plot(T,YNplus(zone,:),'m')
        plot(T,YNmoins(zone,:),'m')
        plot(T,Y(zone,:),'b')
        title(Name_room(zone))
        xlabel("Time(h)")
        ylabel("Temp (°C)")
    end
    legend("Reference","Nominal model","Lower bound","Upper bound","Real model")
end

%% plots 2 Control
if controls
    u_name = ["Bathroom 1";"Bathroom 2";"Bedroom 1";"Bedroom 2";"Bedroom 3";"Bedroom 4";"Lounge"];
    figure('Name',"Control")
    grid on
    for input = 1:7
        subplot(4,2,input)
        plot(T,dU(input,:)*MPCobj.R(input,input)*250,'b')
        hold on
        grid on
        plot(T,U_mpc(input,:)*MPCobj.R(input,input)*250,'r')
        
        title(u_name{input,1})
        xlabel('Temps (h)')
        ylabel('Heat poxer (W)')
    end
    legend("Total input","MPC part of the control")
end

%% Plot 3 condition extérieurs
if scenar_plot
    figure('Name',"Disturbance")
    subplot(3,1,1)
    plot(T,V(4,:),T,Vmd(4,:))
    legend("Real","Predicted")
    title('Exterior Temperature(°)')

    subplot(3,1,2)
    plot(T,V(5,:),T,Vmd(5,:))
    hold on
    plot(T,V(6,:),T,Vmd(6,:))
    legend("Real East","Predicted East","Real West","Predicted West")
    title('Solar Radiation (W/m²)')


    subplot(3,1,3)
    plot(T,V_IG_bath)
    hold on
    plot(T,V_IG_bedroom)
    plot(T,V_IG_kitchen)
    legend("Bathrooms","Bedrooms","Lounge")
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