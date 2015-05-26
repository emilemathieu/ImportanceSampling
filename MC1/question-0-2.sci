// Preliminaries, question 2

function [X]=gauss(N)
  U=rand(1,N,'unif');
  V=rand(1,N,'unif');
  X=sqrt(-2*log(U)) .* sin(2*%pi*V);
endfunction

x=[-4:0.1:4];
densite=(1/sqrt(2*%pi))*exp(- x .* x /2);
plot2d(x,densite,3);

N=10000;
X=gauss(N);
histplot(100,X);

