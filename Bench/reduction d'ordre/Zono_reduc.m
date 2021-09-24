addpath(genpath('..'))

load("full_li_model_coloc",'syslin')
ordre = 12;

[sysb,g,Tt,Ti] = balreal(syslin);
%disp(sum(g(ordre+1:end))*2)
drsys_tc = modred(sysb,ordre+1:size(syslin.A,1),'Truncate');


Tt_1 = Tt(1:ordre,:);
Tt_2 = Tt(ordre+1:end,:);

Ti_1 = Ti(:,1:ordre);
Ti_2 = Ti(:,ordre+1:end);


E = Tt_1*syslin.A*Ti_2*Tt_2;
F = syslin.C*Ti_2*Tt_2;



%pertub set
Zw.c = mean(X2',2);
Zw.R = diag(max(X2',[],2)-min(X2',[],2));
Z = zonotope([Zw.c,Zw.R]);

% Xz = X2;
% Z = ellipsoid.enclosePoints(Xz');
% Z = zonotope(Z);
% Zw.c = Z.Z(:,1);
% Zw.R = Z.Z(:,2:end);
% 
% Xz = X2(:,66:133);
% Z = zonotope.enclosePoints(Xz');
% Zw.c = [Zw.c;Z.Z(:,1)];
% Zw.R = blkdiag(Zw.R,Z.Z(:,2:end));


%initial set
x0 = X2(end,:);
x0 = x0';
Z0.c = Tt_1*x0;
Z0.R = zeros(ordre);

%precision
n = 20;

Zk.c = Z0.c;
Zk.R = Z0.R;
Ts = 15*60; %discretization toute les 15 minutes
t = 0:Ts/3600:72;
for k=1:length(t)
    %méthode 1
    Zk1.c = Zk.c;
    Zk1.R = Zk.R;
    Zk.c = drsys_tc.A*Zk1.c + drsys_tc.B*[Uech(:,k);V2(:,k)] + E*Zw.c;
    Zk.R = [drsys_tc.A*Zk1.R, E*Zw.R];
    Zk.R = reduction(Zk.R, n);
    Zy.c = drsys_tc.C*Zk.c + F*Zw.c;
    Zy.R = [drsys_tc.C*Zk.R, F*Zw.R];
end


    figure; drawZ2D(Zy.c,Zy.R,'m');
    xlabel("temperature zone 1")
    ylabel("temperature zone 2")
%     title(num2str(k)); grid;
    hold off; pause(1);
