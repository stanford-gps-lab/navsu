function p = polyfit2(x,y,n)

% Same as MATLAB polyfit but without some features and the super slow error
% checking. 
x = x(:);
y = y(:);

if nargout > 2
   mu = [mean(x); std(x)];
   x = (x - mu(1))/mu(2);
end

% Construct Vandermonde matrix.
V(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
   V(:,j) = x.*V(:,j+1);
end

% Solve least squares problem.
[Q,R] = qr(V,0);
p = R\(Q'*y);    % Same as p = V\y;


p = p.';          % Polynomial coefficients are row vectors by convention.
