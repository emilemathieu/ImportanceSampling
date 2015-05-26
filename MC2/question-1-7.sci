// Question 5
T=1; // an
S_0=100; 
r=0.05; // par an
sigma=0.3; // par racine d'annee sigma^2 * T  est sans dimension
N=10000

function []=test_call_girsanov(lambda,N)
  W_T = sqrt(T)*rand(1,N,"gauss");
  S_T = S_0*exp((r-sigma^2/2)*T + sigma*(W_T+lambda* T));
  payoff = exp(-r*T) * call(S_T,K);
  importance = exp(-lambda*W_T-lambda**2*T/2) ;
  payoff = importance .* payoff; 
  
  estimation=mean(payoff);  // estimation de la moyenne
  ecart_type=stdev(payoff); // estimation de l'ecart type
  erreur=1.96*ecart_type/sqrt(N); // demi-largeur de l'intervalle de confiance

  printf("Girsanov, lambda=%f, N=%d, %f +- %f\n",lambda, N, estimation, erreur);
endfunction

K=150;
lambda= (log(K/S_0)-(r-sigma^2/2)*T)/(sigma*T)
test_call_girsanov(0,N)
test_call_girsanov(lambda,N)

K=200;
lambda= (log(K/S_0)-(r-sigma^2/2)*T)/(sigma*T)
test_call_girsanov(0,N)
test_call_girsanov(lambda,N)

function [y]=N(x)
  [y,Q]=cdfnor("PQ",x,0,1);
endfunction

function [res] = BS_Call(S_0,K,sigma,r,T)
  d1=(log(S_0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T));
  d2=d1-sigma*sqrt(T);
  res=S_0*N(d1)-K*exp(-r*T)*N(d2);
endfunction
