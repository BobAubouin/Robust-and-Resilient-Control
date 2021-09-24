MPCobj.iMV = 1;
MPCobj.iMO = 1:2;
MPCobj.iUD = [];
MPCobj.iMD = [];

MPCobj.Plant = ss(A,B,eye(2),[],1);
MPCobj.PredictionHorizon = N; %

MPCobj.R = R;
MPCobj.dR = 0;
MPCobj.Q = diag([1,1]);
MPCobj.P = P;

zono.cu = Zu.Z(:,1);
zono.Ru = Zu.Z(:,2:end);

Xmax = 5;
Xmin = -5;
for i =1:2
    MPCobj.OutputVariables(i).Min = Xmin;
    MPCobj.OutputVariables(i).Max = Xmax;
end

zono.cx = Zx.Z(:,1);
zono.Rx = Zx.Z(:,2:end);

zono.cterm = Zterm.Z(:,1);
zono.Rterm = Zterm.Z(:,2:end);
