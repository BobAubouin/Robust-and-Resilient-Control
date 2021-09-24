%% Deny of Service
% Durant cette attaque les données ne sont pas disponible pour une proba
% pu (respectivement py pour les capteurs) à chaque instant, si elle ne sont pas dispo elles sont remplacé par
% les dernière données dispo.

% paramètres de l'attaque DoS
Data.DoS.Tstart = 35; % temps de début de l'attaque (s)
Data.DoS.Tend = 40; % temps de fin de l'attaque (s)

Dos = zeros(length(T),1);
Dos(20:28) = ones(9,1);

Dos(101:108) = ones(8,1);

Data.DoS.du = 0; %strategie de commande de l'actionneur : 0 ne fait rien s'il ne recois rien, 1 reste sur la dernière commande reçu
Data.DoS.dy = 0; %strategie de mesure du controleur : 0 ne fait rien s'il ne recois rien, 1 reste sur la dernière commande reçu

%% Bias injection
% Pour cette attaque, on ajoute progressivement un biais sur la commande ou
% sur les mesures pendant un certain temps.

Data.BI.Tstart = 24*3600; % temps de début de l'attaque (s)
Data.BI.Tend = 36*3600; % temps de fin de l'attaque (s)

Data.BI.BiasU0 = zeros(size_input,1); % Début du biais à 0 pour le modifier doucement et être furtif
Data.BI.transientU = 0.9; %paramètre réglant la vitesse du biai (1=infini, 0=instantané)
Data.BI.FinalBiasU = ones(size_input,1); %valeur finale du biais

Data.BI.BiasY0 = zeros(size_outputs,1); % Début du biais à 0 pour le modifier doucement et être furtif
Data.BI.transientY = 0.9; %paramètre réglant la vitesse du biai (1=infini, 0=instantané)
Data.BI.FinalBiasY = -10*ones(size_outputs,1); %valeur finale du biais


%% Upper Saturation
%Pour cettte attaque on fait saturer tout les actionneurs en espérant tout
%casser

Data.US.Tstart = 45; % temps de début de l'attaque (s)
Data.US.Tend = 55; % temps de fin de l'attaque (s)

Data.US.U = [0;100]; %Max des actionneurs

%% Replay attack
% ici on rejoue un jeu de données précédemment enregistrer pendant que l'on
% réalise une attaque par saturation :

Data.RA.Tstart_save = 12*3600; %Début d'enregistrement (s)
Data.RA.Tend_save = 24*3600; %fin d'enregistrement (s)

Data.RA.Y_save0 = zeros(size_outputs,(Data.RA.Tend_save - Data.RA.Tstart_save)/Ts+1);

Data.RA.Tstart_atk = 36*3600; % temps de début de l'attaque (s)
Data.RA.Tend_atk = 48*3600; % temps de fin de l'attaque (s)

Data.RA.U = [1;1;150;150;150]; %Attaque par saturations

%% False disturbance prediction 

Data.FDP.Tstart_save = 50; %Début d'enregistrement (s)
Data.FDP.Tend_save = 60; %fin d'enregistrement (s)

Data.FDP.Bias = [-5;0]; %Attaque par saturations


%% variable
variable0.RAYsave = Data.RA.Y_save0;
variable0.BIU = Data.BI.BiasU0;
variable0.BIY = Data.BI.BiasY0;

