//**********************************************//
//**************Projet*Math*Fi******************//
//**********************************************//
funcprot(0)

// Pour simuler W a l'instant T on ecrit une fonction
function [W_T] = brownien(N,T,Sigma)
  n=size(Sigma);n=n(2);  
  W_T =  sqrt(T) * Sigma * rand(n,N,"gauss");
endfunction

function [S_T] = black_scholes(N,T,S_0, r,sigma, Sigma)
  n=size(Sigma);n=n(2);
  d=prod(size(S_0));
  S_T = diag(S_0) * (exp(T * diag(r - sigma.^2/2) * ones(d,N) + brownien(N,T,Sigma)));
endfunction

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

function []=test_simple(N,K,fonction_payoff)
  S_T=black_scholes(N,T,S_0, r,sigma, Sigma);
  payoff=exp(-r*T) * fonction_payoff(a,S_T,K);
  estimation=mean(payoff);
  ecart_type=stdev(payoff); // estimation de l'ecart type
  erreur=1.96*ecart_type/sqrt(N); // demi-largeur de l'intervalle de confiance
  printf("Direct     N=%d, %f +- %f\n",N, estimation, erreur);
endfunction

function []=test_avec_vc(N,K,fonction_payoff)
  S_T=black_scholes(N,T,S_0, r,sigma, Sigma);
  I_T=a' * S_T;
  payoff= exp(-r*T) * (fonction_payoff(a,S_T,K)-I_T); // A COMPLETER -- comment utiliser la variable de controle ?
  estimation=mean(payoff) + I_0;// l'esperance de la variable de controle vaut I_0
  ecart_type=stdev(payoff); // estimation de l'ecart type
  erreur=1.96*ecart_type/sqrt(N); // demi-largeur de l'intervalle de confiance
  printf("Direct VC  N=%d, %f +- %f\n",N, estimation, erreur);
endfunction

function []=test_importance_sampling(N,K,fonction_payoff,m)
  W_T = sqrt(T) * rand(d,N,"gauss");
  S_T = diag(S_0) * (exp(T * diag(r - sigma.^2/2) * ones(d,N) + (Sigma * (W_T + m) ) ) );
  payoff = exp(-r*T) * fonction_payoff(a,S_T,K);
  importance = exp(-sum(W_T.*m,1)-sum(m.*m,1)*T/2) ;
  payoff = importance .* payoff; 
  estimation=mean(payoff);  // estimation de la moyenne
  ecart_type=stdev(payoff); // estimation de l'ecart type
  erreur=1.96*ecart_type/sqrt(N); // demi-largeur de l'intervalle de confiance
  printf("Importance N=%d, %f +- %f, mean(drift)=%f\n", N, estimation, erreur, mean(m));
endfunction

function [m]=Find_optimal_drift(K,fonction_payoff)
    m_old=zeros(d,1);
    for i=1:100
        //disp(i)
        W_T = sqrt(T) * Sigma * m_old;
        S_T = diag(S_0) * (exp(T * (r - sigma.^2/2) + W_T ) );
        //disp(fonction_payoff(a,S_T,K)>0)
        for j=1:d
            m(j)=(fonction_payoff(a,S_T,K)>0)*sqrt(T)*(a'*(S_T.*Sigma(:,j)))/(a'*S_T);
        end
        if norm(m-m_old)<0.0005
            break;
        end
        m_old=m;
    end
    m=repmat(m,1,N);
endfunction

//////////// Main Function //////////

//// Initialisation/Calcul des coefficients

// dimension du vecteur d'actifs risqués
d=10;

// Matrice de covariance des mouvements browniens
rho=0.5;
Rho = (1 - rho) *eye(d,d) + rho * ones(d,d); // Matrice de corrélationd es browniens

// Les volatilites de chaque actif
sigma=0.3*ones(d,1);
S_0=100*ones(d,1);

Gamma = diag(sigma) * Rho * diag(sigma); // fabrique la matrice sigma_i * Rho_ij * sigma_j (et oui!)
Sigma=sqroot(Gamma); // Calculons une racine carree de Gamma
//disp(norm(Sigma*Sigma' - Gamma))

T=1; //1 an
r=0.05; //Taux instantanné annuel

// Pondération du panier/de l'indice
a=(1/d)*ones(d,1);

S_0=100*ones(d,1);
I_0= a'*S_0;

//// Appel des fonctions

stacksize(268435454);
N=10000; //Nombre de termes pour la méthode de Monte Carlo


K= 1.5 * a'*S_0; // option out of the money
m=Find_optimal_drift(K,payoff_call_basket);
test_simple(N,K,payoff_call_basket);
test_avec_vc(N,K,payoff_call_basket); // la variable de controle augmente la variance !
test_importance_sampling(N,K,payoff_call_basket,m);

K= 1.0 * a'*S_0; // option at the money
m=Find_optimal_drift(K,payoff_call_basket);
test_simple(N,K,payoff_call_basket);
test_avec_vc(N,K,payoff_call_basket); // la variable de controle marche mais pas terrible
test_importance_sampling(N,K,payoff_call_basket,m);

K= 0.8 * a'*S_0; // option in the money
m=Find_optimal_drift(K,payoff_call_basket);
test_simple(N,K,payoff_call_basket);
test_avec_vc(N,K,payoff_call_basket); 
test_importance_sampling(N,K,payoff_call_basket,m);
    // ce coup ci la variable  de controle sert a qq chose !

// Plus K est petit et mieux ca marche
K= 0.5 * a'*S_0;
m=Find_optimal_drift(K,payoff_call_basket);
test_simple(N,K,payoff_call_basket);
test_avec_vc(N,K,payoff_call_basket); 
test_importance_sampling(N,K,payoff_call_basket,m);