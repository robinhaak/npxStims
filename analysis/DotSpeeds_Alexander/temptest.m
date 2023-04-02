figure
hold on
lambda = 10;
k = 0:100;
%ppoisson = lambda.^k*exp(-lambda)./factorial(k);
lpoiss = -log(poisscdf(k,lambda,'upper'));
lgauss = -log(normcdf(k,lambda,sqrt(lambda),'upper'));
plot(k,lpoiss,'k-')
%plot(k,lgauss,'r-')

%approx =  -exp(-lambda) * (lambda.^(k+1)) ./ sqrt(2*pi*(k+1)).* (k+1/exp(1)).^(k+1)
approx =  (-k+lambda) - k .*log(lambda./k) - log(k+1) + ...
    log(sqrt(2*pi*k)) + log(k+1-lambda);    
plot(k,approx,'r-')


