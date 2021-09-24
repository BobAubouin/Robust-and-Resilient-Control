function [u] = control(MPCobj,xmpc,dy,ref,Vm)
   u = mpcmove(MPCobj,xmpc,dy,ref,Vm);
   q1=0.2; %discrétisation de la commande des volets
   u(1:2)= q1 * round(u(1:2)/q1);
   q2=5; %discrétisation de la commande des radiateurs
   u(3:5)= q2 * round(u(3:5)/q2);
end

