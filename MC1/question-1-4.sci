// Question 4

function [res]=relative_precision(K,N)
  W_T=sqrt(T)*rand(1,N,"gauss");
  S_T=S_0*exp((r-sigma^2/2)*T + sigma*W_T);
  payoff=exp(-r*T) * call(S_T,K);

  estimation=mean(payoff);  // estimation for the price
  ecart_type=st_deviation(payoff); // estimation for the standard deviation
  erreur=1.96*ecart_type/sqrt(N); // half length of the confidence intervall

  res=erreur/estimation;
endfunction

S_0=100;N=10000;
relative_precision(100,N)
relative_precision(150,N)
relative_precision(200,N)
relative_precision(250,N)

