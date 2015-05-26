// Question 2
function [y]=call(x,K)
  y=max(x-K*ones(x),0);
endfunction


function []=test_call(N)
  W_T=sqrt(T)*rand(1,N,"gauss");
  S_T=S_0*exp((r-sigma^2/2)*T + sigma*W_T);
  payoff=exp(-r*T) * call(S_T,K);

  estimation=mean(payoff);  // estimation for the price
  ecart_type=st_deviation(payoff); // estimation for the standard deviation
  method_error=1.96*ecart_type/sqrt(N); // half-width of the confidence interval

  printf("Direct computation N=%d, %f +- %f\n",N, estimation, method_error);
endfunction

test_call(100);
test_call(1000);
test_call(10000);
test_call(100000);

