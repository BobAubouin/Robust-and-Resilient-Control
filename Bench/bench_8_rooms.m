%init
clc;
close all;
clear all;

b_plot = 0; % control the plot of the figures

BRCMRootPath = 'C:\Users\Bob\Documents\IMS\Scripts\tbxmanager\toolboxes\brcm\v1.03\all\brcm\v1.03\all\BRCMToolbox_v1.03';
addpath(genpath(BRCMRootPath));


benchmark = 'benchmark';

thermalModelDataDir =  [pwd,'.\Data_bench_coloc\ThermalModel'];
EHFModelDataDir =      [pwd,'.\Data_bench_coloc\EHFM'];

%% --------------------------------------------------------------------------------------
% 1) Create a building
% --------------------------------------------------------------------------------------

% Create an empty Building object with an optional identifier argument.
buildingIdentifier = benchmark;
B = Building(buildingIdentifier);

%% --------------------------------------------------------------------------------------
% 2) Load the thermal model data
% --------------------------------------------------------------------------------------

% Load the thermal model data. 
B.loadThermalModelData(thermalModelDataDir);

% The thermal model data consists of zones, building elements, constructions,
% materials, windows and parameters. The data of each element group must
% be provided by a separate .xls files and all base files are required for 
% loading the builing data. We require the file names and the file contents to follow a
% specific convention, see the Documentation.

%% --------------------------------------------------------------------------------------
% 3) Declare external heat flux models that should be included
% --------------------------------------------------------------------------------------

% Heat exchange with ambient air and solar gains
EHFModelClassFile = 'BuildingHull.m';                                         % This is the m-file defining this EHF model's class.
EHFModelDataFile = [EHFModelDataDir,filesep,'buildinghull'];                  % This is the spreadsheet containing this EHF model's specification.
EHFModelIdentifier = 'BuildingHull';                                          % This string identifies the EHF model uniquely
B.declareEHFModel(EHFModelClassFile,EHFModelDataFile,EHFModelIdentifier);

% InternalGains
EHFModelClassFile = 'InternalGains.m'; 
EHFModelDataFile = [EHFModelDataDir,filesep,'internalgains']; 
EHFModelIdentifier = 'IG';
B.declareEHFModel(EHFModelClassFile,EHFModelDataFile,EHFModelIdentifier);


% Radiators
EHFModelClassFile = 'Radiators.m'; 
EHFModelDataFile = [EHFModelDataDir,filesep,'radiators']; 
EHFModelIdentifier = 'Rad';
B.declareEHFModel(EHFModelClassFile,EHFModelDataFile,EHFModelIdentifier);

%% --------------------------------------------------------------------------------------
% 4) Display thermal model data to Command Window and draw Building (optional) 
% --------------------------------------------------------------------------------------
if b_plot
    % Print the thermal model data in the Command Window for an overview
    B.printThermalModelData;

    % 3-D plot of Building
    B.drawBuilding;

    % It is possible to control the labeling
    % B.drawBuilding('NoBELabels');
    % B.drawBuilding('NoLabels');
    % B.drawBuilding('NoZoneLabels');

    % 2-D plot of Building
    B.drawBuilding('Floorplan');

    % Drawing parts of the building can be done by a cell array of zone group and/or zone identifiers
    %B.drawBuilding({'ZoneGrp_West'}); 
    B.drawBuilding({'Z0001'}); 
end

%% --------------------------------------------------------------------------------------
% 5) Generate thermal model and full model
% --------------------------------------------------------------------------------------
B.loadThermalModelData(thermalModelDataDir); 
% Generate thermal model (optional)
B.generateThermalModel;

% Generate (full) building model (includes thermal model generation if not yet done)
B.generateBuildingModel;

% Display all available identifiers (these are the names of the control inputs / disturbances / states in the same order as they appear in the matrices)
B.building_model.printIdentifiers;

% Disretization
Ts_hrs = 0.1;
B.building_model.setDiscretizationStep(Ts_hrs);
B.building_model.discretize();

%% --------------------------------------------------------------------------------------
% 7) Retrieve Matrices and generate costs and constraints
% --------------------------------------------------------------------------------------

% Access of full model matrices
discreteTimeFullModelMatrix_A = B.building_model.discrete_time_model.A; % same for Bu,Bv,Bvu,Bxu
continuousTimeFullModelMatrix_A = B.building_model.continuous_time_model.A; % same for Bu,Bv,Bvu,Bxu

% Access of thermal model matrices
continuousTimeThermalModelMatrix_A = B.building_model.thermal_submodel.A; % same for Bq
[discreteTimeThermalModelMatrix_A,discreteTimeThermalModelMatrix_Bq] = B.building_model.thermal_submodel.discretize(Ts_hrs); % these are usually not used, hence not stored

% Access of EHF model matrices (here for the first EHF model)
EHFM1Matrix_Aq = B.building_model.EHF_submodels{1}.Aq; % same for Bq_u, Bq_v, Bq_vu, Bq_xu


% Get constraint matrices such that % Fx*x+Fu*u+Fv*v <= g. These are the constraints for one particular set of potentially 
% time-varying constraintsParameters. Every row of the matrices represents one constraint the name of which is the 
% corresponding entry in constraint_identifiers. The parameters that have to be passed must be in the form 
% constraintsParameters.<EHF_identifier>.<parameters>. Check the documentation to learn which <parameters> are necessary for
% a particular EHF model.


constraintsParameters = struct();


constraintsParameters.BuildingHull.BPos_blinds_E_min = 0.1;
constraintsParameters.BuildingHull.BPos_blinds_E_max = 1;
constraintsParameters.BuildingHull.BPos_blinds_W_min = 0.1;
constraintsParameters.BuildingHull.BPos_blinds_W_max = 1;


constraintsParameters.Rad.Q_rad_Bedroom1_min = 0;
constraintsParameters.Rad.Q_rad_Bedroom1_max = 3;
constraintsParameters.Rad.Q_rad_Bedroom2_min = 0;
constraintsParameters.Rad.Q_rad_Bedroom2_max = 3;
constraintsParameters.Rad.Q_rad_Bedroom3_min = 0;
constraintsParameters.Rad.Q_rad_Bedroom3_max = 3;
constraintsParameters.Rad.Q_rad_Bedroom4_min = 0;
constraintsParameters.Rad.Q_rad_Bedroom4_max = 3;

constraintsParameters.Rad.Q_rad_Lounge_min = 0;
constraintsParameters.Rad.Q_rad_Lounge_max = 3.5;
constraintsParameters.Rad.Q_rad_Bathroom_min = 0;
constraintsParameters.Rad.Q_rad_Bathroom_max = 4;
constraintsParameters.Rad.Q_rad_Bathroom2_min = 0;
constraintsParameters.Rad.Q_rad_Bathroom2_max = 3;

[Fx,Fu,Fv,g] = B.building_model.getConstraintsMatrices(constraintsParameters);


% Get cost vector such that J = cu*u. This is the cost for one particular set of potentially 
% time-varying costParameters. The parameters that have to be passed must be in the form 
% costParameters.<EHF_identifier>.<parameters>. Check the documentation to learn which <parameters> are necessary for
% a particular EHF model. 

costParameters = struct();
costParameters.Rad.costPerJouleHeated = 10;
cu = B.building_model.getCostVector(costParameters);

% If the building_model B.building_model should be saved to use the model in another place, it is necessary that the Classes folder 
% is on the path, otherwise the saved data can not be loaded correctly. If only the matrices are needed, then just 
% the B.building_model.discrete_time_model should be saved and the Classes folder is not necessary.

%% save model
save("discrete_model_coloc_3",'B')


