function [lnL,dL_dTh] = gpr_likelihood(X,Y,theta)

%-- Parameters of the GPR
[n,N] = size(X);
scale = theta(3:end);                                                       % Length-scale
sigmaF2 = theta(2);                                                         % Kernel variance
sigmaU2 = theta(1);                                                         % Noise variance

%-- Constructing the kernel of the training set
K = sigmaF2*sqexp_kern( X, X, scale );
Ky = K + sigmaU2*eye(N);

%-- Cholesky decomposition of the kernel matrix
L = chol(Ky,'lower');
alpha = L'\(L\Y(:));

%-- Calculate the marginal likelihood
lnL = 0.5*( Y(:)'*alpha ) + trace( log(L) );

if nargout == 2
   
    %-- Matrix of partial derivatives of the kernel wrt the hyperparameters
    D = zeros(N,N,n+2);
    D(:,:,1) = eye(N);
    D(:,:,2) = K/sigmaF2;
    for i=1:n
        D(:,:,i+2) = -0.5*( X(i,:) - X(i,:)' ).^2.*K;
    end
    
    dL_dTh = zeros(n+2,1);
    for i=1:n+2
        M1 = (alpha*alpha')*D(:,:,i);
        M2 = L'\(L\D(:,:,i));
        dL_dTh(i) = -0.5*trace( M1 - M2 );
    end
    
end