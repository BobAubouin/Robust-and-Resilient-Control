function J = cost_calcul(Y,U,mpc,mpcdata)
Q = mpc.Q;
R = mpc.R;
dR = mpc.dR;

J = [trace((Y-mpc.Plant.C*mpcdata.x0)'*Q*(Y-mpc.Plant.C*mpcdata.x0)) , trace((U-mpcdata.u0)'*R*(U-mpcdata.u0)) , trace(diff(U')*dR*diff(U')')];
end