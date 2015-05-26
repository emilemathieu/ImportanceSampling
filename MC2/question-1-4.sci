
function [Y] = payoff_call_basket(a,S_T,K)
  // S_T est une matrice d*N
  // I_T et Y sont des vecteurs 1*N
  I_T= a' * S_T;
  Y=max(I_T-K,0);
endfunction


function [Y] = payoff_put_basket(a,S_T,K)
  // S_T est une matrice d*N
  // I_T et Y sont des vecteurs 1*N
  I_T= a' * S_T;
  Y=max(K-I_T,0);
endfunction

function []=test_call_sans_vc(N,K)
  r=0.05;
  T = 1;
  d=10;
  rho=0.5;
  sigma=0.3*ones(d,1);
  S_0=100*ones(d,1);
  Rho = (1 - rho) *eye(d,d) + rho * ones(d,d);
  Gamma = diag(sigma) * Rho * diag(sigma);
  Sigma=sqroot(Gamma);
  a=(1/d)*ones(d,1);
  I_0= a'*S_0;
  
  S_T=black_scholes(N,T,S_0, r,sigma, Sigma);
  payoff=exp(-r*T) * payoff_call_basket(a,S_T,K);
  //disp(payoff)
  estimation=mean(payoff);
  ecart_type=stdev(payoff); // estimation de l'ecart type
  erreur=1.96*ecart_type/sqrt(N); // demi-largeur de l'intervalle de confiance

  printf("Direct N=%d, %f +- %f\n",N, estimation, erreur);
endfunction

function []=test_call_avec_vc(N,K)
  r=0.05;
  T = 1;
  d=10;
  rho=0.5;
  sigma=0.3*ones(d,1);
  S_0=100*ones(d,1);
  Rho = (1 - rho) *eye(d,d) + rho * ones(d,d);
  Gamma = diag(sigma) * Rho * diag(sigma);
  Sigma=sqroot(Gamma);
  a=(1/d)*ones(d,1);
  I_0= a'*S_0;
  
  S_T=black_scholes(N,T,S_0, r,sigma, Sigma);
  I_T=a' * S_T;
  payoff= exp(-r*T) * (payoff_call_basket(a,S_T,K)-I_T); // A COMPLETER -- comment utiliser la variable de controle ?
  estimation=mean(payoff) + I_0;// l'esperance de la variable de controle
                                // vaut I_0
  ecart_type=stdev(payoff); // estimation de l'ecart type
  erreur=1.96*ecart_type/sqrt(N); // demi-largeur de l'intervalle de confiance

  printf("Direct N=%d, %f +- %f\n",N, estimation, erreur);
endfunction

stacksize(1000000);

d=10;
a=(1/d)*ones(d,1);
S_0=100*ones(d,1);

K= 1.0 * a'*S_0; // option at the money
test_call_sans_vc(10000,K);
test_call_avec_vc(10000,K); // la variable de controle marche mais pas terrible

K= 1.2 * a'*S_0; // option out of the money
test_call_sans_vc(10000,K);
test_call_avec_vc(10000,K); // la variable de controle augmente la variance !

K= 0.8 * a'*S_0; // option in the money
test_call_sans_vc(10000,K);
test_call_avec_vc(10000,K); 
    // ce coup ci la variable  de controle sert a qq chose !

// Plus K est petit et mieux ca marche
K= 0.5 * a'*S_0;
test_call_sans_vc(10000,K);
test_call_avec_vc(10000,K); 


