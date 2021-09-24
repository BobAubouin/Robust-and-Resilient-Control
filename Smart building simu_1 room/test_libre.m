n = 7;
x0 = ones(n,1)*0;
x0(1) = -2500;
N = 100;
T = 1:N;

X = zeros(n,N);
X(:,1) = x0;
for i = 1:N-1
    X(:,i+1) = A*X(:,i)+Bu*[1;150];
end

for i=1:n
    subplot(2,4,i)
    plot(T,X(i,:))
end
