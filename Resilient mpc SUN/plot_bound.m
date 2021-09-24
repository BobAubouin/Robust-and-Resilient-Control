% This script plot the norm of the matrix and here bound it is usefull to check wether the paramters h and lambda are optimal or not
% 
% all the paramters needs to be defined first

%% plot bound

T=chol(P);
L_AT = zeros(N+1,1);
L_ABKT = zeros(N+1,1);
B_A = zeros(N+1,1);
B_ABK = zeros(N+1,1);
iT = inv(T);

for i=0:10
    L_AT(i+1) = norm((T*A*iT)^i,2);
    L_ABKT(i+1) = norm((T*(A+B*K)*iT)^i,2);
     B_A(i+1) = h1*lambda1^i;
    B_ABK(i+1) = h2*lambda2^i;
end

figure()
plot(0:N,L_AT,'b');
hold on
plot(0:N,B_A,'r');
xlabel('k',"Interpreter","Latex")
legend('$||TA^kT^{-1}||_2$','$h_1 \lambda_1^k -$',"Interpreter","Latex")
grid on

figure()
plot(0:N,L_ABKT,'b');
hold on
plot(0:N,B_ABK,'r');
xlabel('k',"Interpreter","Latex")
legend('$||T(A+BK)^kT^{-1}||_2$','$h_2 \lambda_2^k -$',"Interpreter","Latex")
grid on
% h1_tilde = h1;% * norm(T,2)*norm(inv(T),2);
% h2_tilde = h2;% * norm(T,2)*norm(inv(T),2);
% 
% Mmax = (-N*log(lambda2) - log(h2))/(log(h1*h2*lambda1/lambda2))