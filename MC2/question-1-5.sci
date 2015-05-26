function [y]=erf(x)
  [y,Q]=cdfnor("PQ",x,0,1);
endfunction

function [res] = BS_Call(S_0,K,sigma,r,T)
  d1=(log(S_0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T));
  d2=d1-sigma*sqrt(T);
  res=S_0*erf(d1)-K*exp(-r*T)*erf(d2);
endfunction

function [y] = black_scholes(lambda, moyenne,var)
  m = moyenne + log(lambda);
  d = (m-log(K))/sqrt(var);
  y= exp(m + (var/2)) * erf(d+sqrt(var)) - K*erf(d)
endfunction

function [I_T,Z_T] = ...
	  simulation_variable_de_controle(r,sigma,Sigma,S_0,a)
  // Simulation de I_T et Z_T 
  n=size(Sigma);n=n(2);
  d=prod(size(S_0));
  LOG_T = T * diag(r - sigma^2/2) * ones(d,N) ...
			 + sqrt(T) * Sigma * rand(n,N,"gauss");
  S_T = diag(S_0) * exp(LOG_T);
  I_T= a' * S_T;
  H_0 = a .* S_0 / I_0; // = a_i S_0^i / I_0
  Z_T = ...; // A COMPLETER
endfunction

function [] = test_call_basket(d,N)

  r=0.05;
  T = 1;
  rho=0.5;
  sigma=0.3*ones(d,1);
  S_0=100*ones(d,1);
  Rho = (1 - rho) *eye(d,d) + rho * ones(d,d);
  Gamma = diag(sigma) * Rho * diag(sigma);
  Sigma=sqroot(Gamma);

  a=(1/d)*ones(d,1);

  I_0 = a' * S_0;
  K = 1.2 * I_0;

  [I_T,Z_T] = simulation_variable_de_controle(r,sigma,Sigma,S_0,a);
  
  // lambda exp(Z_T) est une approximation de S_T
  // on calcule lambda, la moyenne et la variance de Z_T
  H_0 = a .* S_0 / I_0; // = a_i S_0^i / I_0
  J_0 = H_0 .* sigma;   // = a_i S_0^i sigma_i / I_0 

  lambda = I_0; 
  moyenne= T * (H_0' * (r-sigma^2/2)); // moyenne de Z_T
  var = sum(diag(J_0) * Rho * diag(J_0)); // variance de Z_T
  
  // Methode de Monte Carlo sans reduction de variance
  Y=exp(-r*T) * max(I_T-K,0);
  estimation= mean(Y);
  ecart_type=st_deviation(Y); // estimation de l'ecart type
  erreur=1.96*ecart_type/sqrt(N); // demi-largeur de l'intervalle de confiance
  printf("Sans N=%d, %f +- %f\n",N, estimation, erreur);

  
  // Methode de Monte Carlo avec reduction de variance
  Y=exp(-r*T) * (max(I_T-K,0) - max(lambda * exp(Z_T) - K,0));
  estimation= mean(Y) + exp(-r*T) * black_scholes(lambda,moyenne,var);
  ecart_type=st_deviation(Y); // estimation de l'ecart type
  erreur=1.96*ecart_type/sqrt(N); // demi-largeur de l'intervalle de confiance
  printf("Avec N=%d, %f +- %f\n",N, estimation, erreur);

endfunction 

stacksize(10^8);
d=10;
N=1000;
test_call_basket(d,N)



d=10;
a=(1/d)*ones(d,1);
S_0=100*ones(d,1);
I_0 = a'*S_0; 

K=I_0;// option at the money


