load("discrete_model_solo") % data from BRCM

%full discrete systeme
A_k = B.building_model.discrete_time_model.A; % same for Bu,Bv,Bvu,Bxu
Bu_k = B.building_model.discrete_time_model.Bu;
Bv_k = B.building_model.discrete_time_model.Bv;
Bvu_k = B.building_model.discrete_time_model.Bvu;
Bxu_k = B.building_model.discrete_time_model.Bxu;
C_k = [1,0,0,0,0,0,0];

Ts = 1;
Plant = ss(A_k,[Bu_k,Bv_k],C_k,[],Ts);

%Plant.InputName = cat(1,B.building_model.identifiers.u,B.building_model.identifiers.v);
Plant.OutputName = {'T_zone_1'};

iMV = 1:2;
iUD = 3;
iMD = 4:5;
iMO = 1;

Plant = setmpcsignals(Plant,'MV',iMV,'MD',iMD,'UD',iUD, ...
    'MO',iMO,'UO',[]);

%Ts = 60*15;% 15 minutes

MPCobj = mpc(Plant,Ts);

MPCobj.PredictionHorizon = 8; % 4h
MPCobj.ControlHorizon = MPCobj.PredictionHorizon;

Q = [1000];
R = [10000,1];

MPCobj.Weights.ManipulatedVariables = R;
MPCobj.Weights.OutputVariables = Q;

MPCobj.ManipulatedVariables(1).Min = 0;
MPCobj.ManipulatedVariables(1).Max = 1;
MPCobj.ManipulatedVariables(2).Min = 0;
MPCobj.ManipulatedVariables(2).Max = 100;

Tmax = 40;
Tmin = 0;
for i =1:1
    MPCobj.OutputVariables(i).Min = Tmin;
    MPCobj.OutputVariables(i).Max = Tmax;
end

zono.cx = ones(8,1)*mean([Tmax,Tmin]);
zono.Rx = diag(repmat(mean([Tmax,Tmin])-Tmin,8,1));


setEstimator(MPCobj,'custom');
xmpc = mpcstate(MPCobj);

%% observer
% A = MPCobj.Model.Plant.A;
% C = MPCobj.Model.Plant.C;
% p = -0.89 + (1:size(A,1)/2)/(size(A,1)*20);
% pim = 1i*(1:size(A,1)/2)/(size(A,1)*5);
% p1 = p + pim;
% p2 = p - pim;
% poles = [p1,p2];
% 
% % [L,M,Am,Cm,Bu,Bv,Dvm] = getEstimator(MPCobj);
% % Poles = eig(Am - L*Cm);
% % poles = Poles(1:size(A,1));
% 
% L = place(A',C',poles)';


