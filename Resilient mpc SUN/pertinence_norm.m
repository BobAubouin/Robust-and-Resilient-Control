% Ce script test la pertinence de la borne proposé par SUN dans son article.
% Pour un horizon de prédiction N donné ce script test toute les durées 
% d'attaque pour toute les combinaison possible et  compare la borne de la trajectoire la plus grande avec la borne proposée par SUN

% définition du système
A = [1   -1.2; 1.2 1.1];
B = [1 ; 0.5]; 

K = -[1.5718    0.2611];
P = [1.9385, 1.7088 ; 1.7088, 5.0552];
T = chol(P);
N = 10;

%boucle sur le nombre d'attaque
Normmax=[];
Normmax2=[];
borne = [];
Mmax = 7;
for M=1:Mmax
DoS = create_attack(N,M); %cette fonction renvoie toute les combinaisons d'attaque
[h1,h2,lambda1,lambda2,ok] = get_param_h_l(A,B,K,P,N,M,1);
f = @(M) h1^M*h2^(M+1)*lambda1^M*lambda2^(N-M);
Norm = zeros(size(DoS,2),1);
Norm2 = zeros(size(DoS,2),1);
for k=1:size(DoS,2)
    XN = eye(2);
    XN2 = eye(2);
    for i=1:N
       XN= XN*(A + DoS(i,k)*B*K);   
       XN2= XN2*T*(A + DoS(i,k)*B*K)*inv(T);  
    end
    Norm2(k) = norm(XN2,2);
end
Normmax2=[Normmax2 max(Norm2)];

borne = [borne f(M)];
end

plot(1:Mmax,log(Normmax2))
hold on
plot(1:Mmax,log(borne))
grid on
legend("$ln(||\prod_{k=0}^{N-1} A + \nu_k BK||_P)$","$ln(h_1^M h_2^{M+1}\lambda_1^M \lambda^{N-M})$","Interpreter","Latex",'FontSize',14)
xlabel('M')      