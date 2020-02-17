function H = hermite(x,pb)

x = x(:);
N = length(x);

H = ones(N,pb);
if pb >= 2
    H(:,2) = 2*x;
end
for i=3:pb
    n = i-1;
    H(:,i) = 2*x.*H(:,i-1) - 2*(n-1)*H(:,i-2);
end

for i=1:pb
    n = i-1;
    c = ( 2^n*factorial(n)*sqrt(pi) )^(-1/2);
    H(:,i) = c*H(:,i).*exp(-(x.^2)/2);
end
