function Zset = getRPIset(ABK,D,Zw,alpha)
%This function return an RPI qith the method described in " Invariant Approximations of the Minimal Robust
% Positively Invariant Set" from S. V. Rakovic 2005
% Input :   - ABK : Strictly stable matrix (p(ABK)<1) of the system
%           - D : matrice of the distrubances
%           - Zw: Set of all the possible disturbances (define with Cora
%           toolbox)
%           - alpha: scalar paramete between ]0,1[, the smaller it is, the
%           smaller is the output set
% Output:   - Zset: RPI calculated
%
% Bob Aubouin, 23/06/2021
s = 0;
while ~in(alpha*D*Zw , ABK^s*D*Zw)
    s=s+200;
%         disp(s);
%         plot(alpha*D*Zw);
%         hold on
%         plot(ABK^s*D*Zw);
%         input('continuer ?');
end
%input('continuer ?');
Zset = 0*D*Zw;
for i = 0:s-1
    Zset = Zset + ABK^i*D*Zw;
    if mod(i,2000)==0
        %Zset = reduce(Zset,'combastel',100);
        %disp(i);
    end
end

int = 1/(1-alpha);
Zset = Zset*int;
end
