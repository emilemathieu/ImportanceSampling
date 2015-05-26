// Question 5
function []=test_call_girsanov(S_0,K,sigma,r,T, lambda,NumberOfDrawings)
  W_T = sqrt(T)*rand(1,NumberOfDrawings,"gauss");
  S_T = S_0*exp((r-sigma^2/2)*T + sigma*(W_T+lambda* T));
  payoff = exp(-r*T) * call(S_T,K);
  importance = exp(-lambda * W_T - ((lambda^2*T) /2.0));
  payoff = importance .* payoff;
  
  estimation=mean(payoff); 
  ecart_type=st_deviation(payoff);
  method_error=1.96*ecart_type/sqrt(NumberOfDrawings);

  printf("Girsanov method : number of drawings = %d, %f +- %f\n",  ...
         NumberOfDrawings, estimation, method_error);
endfunction

function [y]=N(x)
  [y,Q]=cdfnor("PQ",x,0,1);
endfunction

function [res] = BS_Call(S_0,K,sigma,r,T)
  d1=(log(S_0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T));
  d2=d1-sigma*sqrt(T);
  res=S_0*N(d1)-K*exp(-r*T)*N(d2);
endfunction

r=0.05;
sigma=0.3;
T=1;

// Works when S_0 < K, especially when S_0 << K
// With K=100 try S_0=50,70,90,100
K=100;

NumberOfDrawings=10000;

// Very efficient in this case
S_0=50; 
lambda= (log(K/S_0)-(r-sigma^2/2)*T)/(sigma*T);
// without importance sampling
test_call_girsanov(S_0,K,sigma,r,T, 0, NumberOfDrawings);
// with importance sampling
test_call_girsanov(S_0,K,sigma,r,T, lambda, NumberOfDrawings);

// Efficient but less than in the first case
S_0=70; 
lambda= (log(K/S_0)-(r-sigma^2/2)*T)/(sigma*T);
// without importance sampling
test_call_girsanov(S_0,K,sigma,r,T, 0, NumberOfDrawings);
// with importance sampling
test_call_girsanov(S_0,K,sigma,r,T, lambda, NumberOfDrawings);

// Even less efficient
S_0=90; 
lambda= (log(K/S_0)-(r-sigma^2/2)*T)/(sigma*T);
// without importance sampling
test_call_girsanov(S_0,K,sigma,r,T, 0, NumberOfDrawings);
// with importance sampling
test_call_girsanov(S_0,K,sigma,r,T, lambda, NumberOfDrawings);


// Nothing happen !
S_0=100; 
lambda= (log(K/S_0)-(r-sigma^2/2)*T)/(sigma*T);
// without importance sampling
test_call_girsanov(S_0,K,sigma,r,T, 0, NumberOfDrawings);
// with importance sampling
test_call_girsanov(S_0,K,sigma,r,T, lambda, NumberOfDrawings);