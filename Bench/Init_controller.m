load("discrete_model_coloc") % data from reduc model
load("reduced_model_coloc",'Tt','drsys_tc_12') % data from reduc model

Plant = drsys_tc_12;

%Plant.InputName = cat(1,B.building_model.identifiers.u,B.building_model.identifiers.v);
Plant.OutputName = {'T_zone_1','T_zone_2','T_zone_3','T_zone_4','T_zone_5','T_zone_6','T_zone_7','T_zone_8'};

iMV = 1:size_input;
iUD = size_input+1:size_input+(size_disturbances-size_predicted_disturbances);
iMD = size_input+(size_disturbances-size_predicted_disturbances) + 1:size_input+size_disturbances;
iMO = 1:size_outputs;

Plant = setmpcsignals(Plant,'MV',iMV,'MD',iMD,'UD',iUD, ...
    'MO',iMO,'UO',[]);

Ts = 60*15;% 15 minutes

MPCobj = mpc(Plant,Ts);

MPCobj.PredictionHorizon = 2; % 4h
MPCobj.ControlHorizon = MPCobj.PredictionHorizon;

Q = [10,10,10,10,0.1,1,10,10]*1000;
P = [0,0,1,1,1,1,1,1,1];
R = [1,1,0.1,0.1,0.1,0.1,0.1,0.1,0.1];

MPCobj.Weights.ManipulatedVariables = P;
MPCobj.Weights.ManipulatedVariablesRate = R;
MPCobj.Weights.OutputVariables = Q;

MPCobj.ManipulatedVariables(1).Min = 0;
MPCobj.ManipulatedVariables(1).Max = 1;
MPCobj.ManipulatedVariables(2).Min = 0;
MPCobj.ManipulatedVariables(2).Max = 1;
MPCobj.ManipulatedVariables(3).Min = 0;
MPCobj.ManipulatedVariables(3).Max = 150;
MPCobj.ManipulatedVariables(4).Min = 0;
MPCobj.ManipulatedVariables(4).Max = 150;
MPCobj.ManipulatedVariables(5).Min = 0;
MPCobj.ManipulatedVariables(5).Max = 150;

zono.cu = [0.5;0.5;75;75;75];
zono.Ru = diag([0.5;0.5;75;75;75]);

Tmax = 30;
Tmin = 10;
for i =1:8
    MPCobj.OutputVariables(i).Min = Tmin;
    MPCobj.OutputVariables(i).Max = Tmax;
end

zono.cx = ones(8,1)*mean([Tmax,Tmin]);
zono.Rx = diag(repmat(mean([Tmax,Tmin])-Tmin,8,1));


setEstimator(MPCobj,'custom');
xmpc = mpcstate(MPCobj);

%% observer
A = MPCobj.Model.Plant.A;
C = MPCobj.Model.Plant.C;
p = -0.89 + (1:size(A,1)/2)/(size(A,1)*20);
pim = 1i*(1:size(A,1)/2)/(size(A,1)*5);
p1 = p + pim;
p2 = p - pim;
poles = [p1,p2];

% [L,M,Am,Cm,Bu,Bv,Dvm] = getEstimator(MPCobj);
% Poles = eig(Am - L*Cm);
% poles = Poles(1:size(A,1));

L = place(A',C',poles)';


