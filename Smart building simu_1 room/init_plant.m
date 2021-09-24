load("discrete_model_solo") % data from BRCM

%full discrete systeme
A_k = B.building_model.discrete_time_model.A; % same for Bu,Bv,Bvu,Bxu
Bu_k = B.building_model.discrete_time_model.Bu;
Bv_k = B.building_model.discrete_time_model.Bv;
Bvu_k = B.building_model.discrete_time_model.Bvu;
Bxu_k = B.building_model.discrete_time_model.Bxu;
C_k = [1,0,0,0,0,0,0];


%% system sizes

size_input = size(Bu_k,2);
size_disturbances = size(Bv_k,2);
size_predicted_disturbances = 2;
size_outputs = size(C_k,1);
size_state = size(Bu_k,1);


