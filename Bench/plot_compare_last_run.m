%% plot erreur de poursuite :
if errors
    figure('Name',"error")
    for zone = 1:8
        subplot(4,2,zone)
        plot(T/3600,Y(zone,:)-ref*ones(1,N),T/3600,Y2(zone,:)-ref*ones(1,N))
        legend("run 1","run 2")
        title(['Ecart de température réelle zone ',num2str(zone)])
        xlabel('Temps (h)')
        ylabel('Température (°c)')
    end
end
%% plots 1 Y
if outputs
    figure('Name',"Outputs")
    for zone = 1:8
        subplot(4,2,zone)
        plot(T/3600,Y(zone,:),'b')
        hold on
        plot(T/3600,Y2(zone,:),'r')
        plot(T/3600,ref*ones(1,N),'g')
        legend("current run","saved run","référence")
        title(['Température zone ',num2str(zone)])
        xlabel('Temps (h)')
        ylabel('Température (°c)')
    end
end
%% plots 2 Command
if controls
    figure('Name',"Control")
    for input = 1:5
        subplot(2,3,input)
        plot(T/3600,dU(input,:),'b')
        hold on
        plot(T/3600,dU2(input,:),'r')
        legend("current run","saved run")
        title(['Commande',num2str(input)])
        xlabel('Temps (h)')
        ylabel('W/m²')
    end
end
%% observer

if observer
    figure('Name',"Observer")
    C = Plant.C;
    Yo = C*Xo;
    Yo2 = C*Xo2;
    for zone = 1:8
        subplot(4,2,zone)
        plot(T/3600,Y(zone,:),'b')
        hold on
        plot(T/3600,Yo(zone,:),'b--')
        plot(T/3600,Y2(zone,:),'r')
        plot(T/3600,Yo2(zone,:),'r--')
        plot(T/3600,ref*ones(1,N),'g')
        legend("real current run","observed current run","real saved run","observed saved run","référence")
        title(['Température zone ',num2str(zone)])
        xlabel('Temps (h)')
        ylabel('Température (°c)')
    end
end

%% zonotope

if zonoplot
   figure('Name',"Zono")
   grid on
   for i =1:N
    drawZ2D(zono.cz+XN(:,i),zono.Rz,'line,b');
    hold on
    drawZ2D(zono2.cz+XN2(:,i),zono2.Rz,'line,g'); 
   end
   legend("Z set last","Z set saved",'Interpreter','latex')
   plot(XN(1,:),XN(2,:),'b*-',"DisplayName","Nominal last")
   plot(XN2(1,:),XN2(2,:),'g*-',"DisplayName","Nominal saved")
   plot(Y(1,:),Y(2,:),'ro-',"DisplayName","Real last")
   plot(Y2(1,:),Y2(2,:),'ko-',"DisplayName","Real saved")
   legend()
   ylabel('$x_2$','Interpreter','latex')
   xlabel('$x_1$','Interpreter','latex')
end