load("discrete_model_coloc") % data from reduc model
load("full_li_model_coloc",'syslin')

Plant = syslin;

Plant.InputName = cat(1,B.building_model.identifiers.u,B.building_model.identifiers.v);
Plant.OutputName = {'T_zone_1','T_zone_2','T_zone_3','T_zone_4','T_zone_5','T_zone_6','T_zone_7','T_zone_8'};

iMV = 1:size_input;
iUD = size_input+1:size_input+(size_disturbances-size_predicted_disturbances);
iMD = size_input+(size_disturbances-size_predicted_disturbances) + 1:size_input+size_disturbances;
iMO = 1:size_outputs;

Plant = setmpcsignals(Plant,'MV',iMV,'MD',iMD,'UD',iUD, ...
    'MO',iMO,'UO',[]);

Ts = 60*15;% 15 minutes

MPCobj = mpc(Plant,Ts);

MPCobj.PredictionHorizon = 4/0.25; % 4h


MPCobj.Weights.ManipulatedVariables = [0,0,1,1,1];
MPCobj.Weights.ManipulatedVariablesRate = [1,1,0.1,0.1,0.1];
MPCobj.Weights.OutputVariables = [10,10,10,10,0,1,10,10];

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

for i =1:8
    MPCobj.OutputVariables(i).Min = 10;
    MPCobj.OutputVariables(i).Max = 30;
end



xmpc = mpcstate(MPCobj);