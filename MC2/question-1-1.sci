d=10;
rho=0.5;

// 1 sur la diagonale
// rho en dehors de la diagonale

Rho = (1 - rho) *eye(d,d) + rho * ones(d,d);

// spec calcule les valeurs propre de Gamma (son spectre)

spec(Rho)

// verifier que toutes les valeurs propres sont positives

// ...

// fabrique la matrice sigma_i * Rho_ij * sigma_j (et oui!)

sigma = 0.3 * ones(1,d); 
Gamma = diag(sigma) * Rho * diag(sigma);

spec(Gamma)

