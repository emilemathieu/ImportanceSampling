// Question 3

// C-P = S_0 - K exp(-rT)
// So we can use a new estimator
//
//    S_0 - K exp(-rT) + exp(-rT) * (K-S_T)_+

function [y]=put(x,K)
  y=max(K*ones(x)-x,0);
endfunction

function []=test_call_arbitrage(N)
  W_T=sqrt(T)*rand(1,N,"gauss");
  S_T=S_0*exp((r-sigma^2/2)*T + sigma*W_T);
  payoff=S_0-K*exp(-r*T) + exp(-r*T) * put(S_T,K);

  estimation=mean(payoff);  // estimation for the price
  ecart_type=st_deviation(payoff); // estimation for the standard deviation
  method_error=1.96*ecart_type/sqrt(N); // half length of the confidence interval

  printf("CallPut duality method N=%d, %f +- %f\n",N, estimation, method_error);
endfunction

K=100;test_call(1000);test_call_arbitrage(1000);
K= 80;test_call(1000);test_call_arbitrage(1000);
K= 60;test_call(1000);test_call_arbitrage(1000);

