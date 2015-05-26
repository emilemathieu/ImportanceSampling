T=1; // a year
S_0=100; 
r=0.05; // by year
sigma=0.3; // by square root of year, sigma^2 * T has no dimension
K=100;

// Question 1
N=100000;
W_T=sqrt(T)*rand(1,N,"gauss");
S_T=S_0*exp((r-sigma^2/2)*T + sigma*W_T);

histplot(100,W_T);
histplot(100,S_T);

