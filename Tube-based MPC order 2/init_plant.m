%% Parameters 
%systeme
Ts = 1;
A = [1   -1.2; 1.2 1.1];
%B = [0.4; 0.2]; %OK
B = [1 ; 0.5]; %NOK

C = eye(2);

%ensembles perturbation taille 1
D = [0.2; 0.4 ];
Zw.c = 0.02;
Zw.R = 0.1;
Zd = zonotope([Zw.c Zw.R]);


%ensembles perturbation taille 2
D = [0.2 0.1; 0.4 0.1];
Zw.c = [0.02; 0];
Zw.R = [0.1 0; 0 0.005];
Zd = zonotope([Zw.c Zw.R]);


%% system sizes

size_input = size(B,2);
size_disturbances = size(D,2);
size_predicted_disturbances = 0;
size_outputs = size(C,1);
size_state = size(A,1);
x0 = [-0.4;0.8];

