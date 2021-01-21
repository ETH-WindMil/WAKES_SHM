function k = sqexp_kern(x,x0,l)

[n,m1] = size(x);
[~,m2] = size(x0);

k = zeros(m1,m2);

if numel(l) == 1
    for i=1:m2
        d = x - repmat(x0(:,i),1,m1);
        r = diag( d'*d );
        k(:,i) = exp( -(r*l^2)/2 );
    end
elseif numel(l) == n
    for i=1:m2
        d = x - repmat(x0(:,i),1,m1);
        r = diag( d'*diag(l.^2)*d );
        k(:,i) = exp( -r/2 );
    end
end