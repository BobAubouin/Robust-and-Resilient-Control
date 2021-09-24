function result = f(M,N,h1,h2,lambda1,lambda2)
if N<=M
    hlambda_m = max(h1*lambda1,h2*lambda2);
    result = hlambda_m^N;
elseif M<=N/2
    result = h1^M*h2^min(M+1,N-M)*lambda1^M*lambda2^(N-M);
else
    result = h1^(N-M+1)*h2^(N-M)*lambda1^M*lambda2^(N-M);
end
end