load("discrete_model_coloc_3") % data from BRCM

%full discrete systeme
A_k = B.building_model.discrete_time_model.A; % same for Bu,Bv,Bvu,Bxu
Bu_k = B.building_model.discrete_time_model.Bu;
Bv_k = B.building_model.discrete_time_model.Bv;
Bvu_k = B.building_model.discrete_time_model.Bvu;
Bxu_k = B.building_model.discrete_time_model.Bxu;
C_k = diag((1:size(A_k,1))<=8);
C_k = double(C_k(1:8,:));


%% system sizes

size_input = size(Bu_k,2);
size_disturbances = size(Bv_k,2);
size_predicted_disturbances = size(Vm.Data,1);
size_outputs = size(C_k,1);
size_state = size(Bu_k,1);
x0 = 18*ones(size_state,1);

