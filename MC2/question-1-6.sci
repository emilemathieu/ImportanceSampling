// Question 4
T=1; // an
S_0=100; 
r=0.05; // par an
sigma=0.3; // par racine d'annee sigma^2 * T  est sans dimension


function [res] = call(x,k)
  res=max(x-k,0); // A COMPLETER
endfunction

function [res]=precision_relative(K,N)
  W_T=sqrt(T)*rand(1,N,"gauss");
  S_T=S_0*exp((r-sigma^2/2)*T + sigma*W_T);
  payoff=exp(-r*T) * call(S_T,K);

  estimation=mean(payoff);  // estimation de la moyenne
  ecart_type=st_deviation(payoff); // estimation de l'ecart type
  erreur=1.96*ecart_type/sqrt(N); // demi-largeur de l'intervalle de confiance

  res=erreur/estimation;
endfunction

S_0=100;N=10000;
precision_relative(100,N)
precision_relative(150,N)
precision_relative(200,N)
precision_relative(250,N)

