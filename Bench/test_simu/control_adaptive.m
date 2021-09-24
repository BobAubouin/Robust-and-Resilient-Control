function [u] = control_adaptive(MPCobj,xmpc,dy,ref,Vm)
   Npred = MPCobj.PredictionHorizon;
   Vm = reshape(Vm,[size(Vm,1),size(Vm,3)]);
   Vpred = Vm(:,1:min(Npred,size(Vm,2)));
   u = mpcmoveAdaptive(MPCobj,xmpc,MPCobj.Model.Plant,[],dy,ref,Vpred');
%    q1=0.2; %discrétisation de la commande des volets
%    u(1:2)= q1 * round(u(1:2)/q1);
%    q2=5; %discrétisation de la commande des radiateurs
%    %u(3:5)= q2 * round(u(3:5)/q2);
end

