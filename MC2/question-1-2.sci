// Pour simuler W a l'instant T on ecrit une fonction
function [W_T] = brownien(N,T,Sigma)

  // Sigma est obtenu par sqroot(Gamma)
  // et lorsque la matrice Gamma est non inversible (ex: rho = 1).
  // et le "sqroot" de scilab renvoie une matrice dxn (n<d) au lieu de dxd
  // ou n est le rang de Gamma
  n=size(Sigma);n=n(2);  

  // apres c'est simple (pourquoi ca marche ?)
  W_T =  sqrt(T) * Sigma * rand(n,N,"gauss");
endfunction

// S_T
function [S_T] = black_scholes(N,T,S_0, r,sigma, Sigma)
  n=size(Sigma);n=n(2);
  d=prod(size(S_0));
  S_T = exp(T * diag(r - sigma.^2/2) * ones(d,N) + brownien(N,T,Sigma));
  S_T = diag(S_0) * S_T;// produit terme a terme
endfunction

d=3;
rho=0.5;
// rho = 1

// Les volatilites de chaque actif

sigma=0.3*ones(d,1);
S_0=100*ones(d,1);

// 1 sur la diagonale
// rho en dehors de la diagonale

Rho = (1 - rho) *eye(d,d) + rho * ones(d,d);

// fabrique la matrice sigma_i * Rho_ij * sigma_j (et oui!)

Gamma = diag(sigma) * Rho * diag(sigma);

// Calculons une racine carree de Gamma
// sqroot fait le travail
// attention il y a d'autre fonction "racine carree" 
// qui ne font pas ce que l'on veut !

Sigma=sqroot(Gamma);

// On verifie on sait jamais, on doit obtenir un truc petit ...

norm(Sigma*Sigma' - Gamma)

N=5;
T=1;
r=0.05;

brownien(N,T,Sigma)'
black_scholes(N,T,S_0, r,sigma, Sigma)'
