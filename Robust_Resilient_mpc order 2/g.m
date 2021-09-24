function result = g(M,N,h1,h2,lambda1,lambda2)
hlambda_m = max(h1*lambda1,h2*lambda2);
if 2*M <= N-3
   result = (hlambda_m-(hlambda_m)^(M+1))/(1-hlambda_m)+(h1*lambda1)^M*(h2^(M+1)*(lambda2^(M+2)-lambda2^(N-M))/(1-lambda2) + ((h2*lambda2)-(h2*lambda2)^(M+2))/(1-h2*lambda2))+1;
elseif 2*M > N-3 && M <= N-2
   result = (hlambda_m-(hlambda_m)^(M+1))/(1-hlambda_m)+(h1*lambda1)^M*((h2*lambda2)-(h2*lambda2)^(N-M))/(1-h2*lambda2)+1;
else
    result = (hlambda_m-(hlambda_m)^(N))/(1-hlambda_m)+1;
end
end
