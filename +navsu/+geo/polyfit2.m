function p = polyfit2(x,y,n)

% Same as MATLAB polyfit but without some features and the super slow error
% checking. 
if isrow(x)
    x = x';
end
if size(y, 1) ~= length(x)
    y = y';
end

if nargout > 2
   mu = [mean(x); std(x)];
   x = (x - mu(1))/mu(2);
end

% Construct Vandermonde matrix.
V = x.^(n:-1:0);

% Solve least squares problem.
[Q, R] = qr(V, 0);
p = R\(Q'*y);    % Same as p = V\y;


p = p.';          % Polynomial coefficients are row vectors by convention.

end