// Preliminaries, question 3

function []=test_1(N,b)
  X=exp(b*rand(1,N,"gauss"));

  estimation=mean(X);
  etimated_error=1.96 * st_deviation(X)/sqrt(N);
  exact_value = exp(b*b/2);
  method_error = abs(estimation-exact_value)/exact_value;
  printf("N=%d, Value = %f Relative error = %f\%", ...
          N,estimation,100 * method_error); 
endfunction

stacksize(10000000);

test_1(100000,1);
test_1(100000,2);
test_1(100000,4);
test_1(100000,6);
test_1(100000,8);
test_1(1000000,8);


