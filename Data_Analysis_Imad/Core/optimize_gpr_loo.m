function [hyperpar,lnL] = optimize_gpr_loo(Xtrain,Ytrain,theta0)

n = size(Xtrain,1);
problem.objective = @(theta) gpr_loo_lnL( Xtrain, Ytrain, theta );
problem.solver = 'fmincon';
problem.lb = zeros(n+2,1);

if nargin == 3
    problem.x0 = theta0;
else
    problem.x0 = 1e-3*ones(1,n+2);
end

problem.options = optimoptions('fmincon');
problem.options.Display = 'iter';
% problem.options.PlotFcns = @optimplotfval;
problem.options.UseParallel = true;
% problem.options.CheckGradients = true;
problem.options.SpecifyObjectiveGradient = true;

[hyperpar,lnL] = fmincon(problem);


%--------------------------------------------------------------------------
function [lnL_loo,dL_dTh] = gpr_loo_lnL( X, Y, theta )

%-- Parameters of the GPR
[n,N] = size(X);
scale = theta(3:end);                                                       % Length-scale
sigmaF2 = theta(2);
sigmaU2 = theta(1);                                                         % Noise variance

%-- Constructing the kernel of the training set
K = sigmaF2*sqexp_kern( X, X, scale );
Ky = K + sigmaU2*eye(N);

%-- Cholesky decomposition of the kernel matrix
alpha = Ky\Y(:);

%-- LOO pseudo-likelihood
J = inv(Ky);

% Calculating the LOO pseudo-likelihood
lnL_loo = 0.5*mean( alpha.^2./diag(J) ) - 0.5*mean( log(diag(J)) ) + 0.5*log(2*pi);

if nargout == 2
    
    %-- Matrix of partial derivatives of the kernel wrt the hyperparameters
    D = zeros(N,N,n+2);
    D(:,:,1) = eye(N);
    D(:,:,2) = K/sigmaF2;
    for i=1:n
        D(:,:,i+2) = -0.5*( X(i,:) - X(i,:)' ).^2.*K;
    end
    
    Z = D;
    M2 = zeros(N,N,n+2);
    for i=1:n+2
        Z(:,:,i) = Ky\D(:,:,i);
        M2(:,:,i) = Z(:,:,i)/Ky;
    end
    
    dL_dTh = zeros(n+2,1);
    
    for j=1:n+2
        sj = diag(M2(:,:,j));
        rj = -M2(:,:,j)*Y(:);
        dL_dTh(j) = 0.5*mean( ( 1 + alpha.^2./diag(J) ).*( sj./diag(J) ) ) ...
            + mean( alpha.*rj./diag(J) );
    end
    
end