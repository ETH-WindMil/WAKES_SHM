function [Yast,sigmaY2] = gpr_predict(X,Y,Xast,theta)

%-- Parameters of the GPR
[~,N] = size(X);
scale = theta(3:end);                                                       % Length-scale
sigmaF2 = theta(2);                                                         % Kernel variance
sigmaU2 = theta(1);                                                         % Noise variance

%-- Constructing the kernel of the training set
K = sigmaF2*sqexp_kern( X, X, scale );
Ky = K + sigmaU2*eye(N);

%-- Cholesky decomposition of the kernel matrix
L = chol(Ky,'lower');
alpha = L'\(L\Y(:));

%-- Kernel evaluated between training and testing points
K_ast = sigmaF2*sqexp_kern( X, Xast, scale );
v = L\K_ast;

%-- Prediction
Yast = K_ast'*alpha;
sigmaY2 = sigmaF2 - diag(v'*v);