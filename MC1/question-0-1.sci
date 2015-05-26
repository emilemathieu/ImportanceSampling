// Preliminaries, question 1

function [res]=Mean(x)
  res=sum(x)/prod(size(x));
endfunction

function [res]=Variance(x)
  res = Mean(x^2)-Mean(x)^2;
  N = prod(size(x));
  res = res*N/(N-1);
endfunction

n=1000;
X=rand(1,n,"unif");
Mean(X)
mean(X)
variance(X)
Variance(X)
