function yi = lagrangeInter(x,y,xi)

% SYNTAX:
%   yi = LagrangeInter(x,y,xi);
%
% INPUT:
%   x,y = row-vectors of (n) data values (x,y)
%       NOTE: y can be an 
%   xi  = a row-vector of x-values, where interpolation is to be found (could be a single value)
%
% OUTPUT:
%   yi = a row-vector of interpolated y-values
%
% DESCRIPTION:
%   Lagrange interpolation algorithm.

%----------------------------------------------------------------------------------------------
%
% Author: Dmitry Pelinovsky
% Available online at: http://dmpeli.mcmaster.ca/Matlab/Math4Q3/Lecture2-1/LagrangeInter.m
% 
% Original code modified by Fabian Rothmaier. Inner loops vectorized for
% improved computational performance. Additional dimension checks added to
% inputs. Vector length now declared as "n" for streamlined notation.
%----------------------------------------------------------------------------------------------

% expect to work with certain vector shapes
if isrow(x)
    x = x';
end
n = length(x); % the degree of interpolation polynomial

if size(y, 2) ~= n
    y = y'; % for dot product at the end
end
if iscolumn(xi)
    xi = xi'; % need this as row vector
end

ni = length(xi); % the number of x-values, where interpolation is to be found

L = ones(n, ni); % the matrix for Lagrange interpolating polynomials L_(n,k)(x)
                    % has (n) rows for each polynomial at k = 1,...,n
                    % has ni column for each x-value of xi
                    
% Note: the algorithm uses the MATLAB capacities for matrix handling!
% The two nested loops below are designed to compute in parallel 
% the values of Lagrange interpolating polynomials at each x-value of xi !
% Edit by Fabian Rothmaier: Matlab capacities leveraged even more. Can
% still be called with a vector of values xi.

for k = 1 : n % start the outer loop through the data values for x
    
    Lk = L(k, :);
    L = L ./ (x - x(k)) .* (xi - x(k));
    L(k, :) = Lk; % reset this one
    
end % the end of the outer loop

% Now multiply the values for Lagrange interpolating polynomials by the data values for y
yi = y * L; % use matrix multiplication of row-vector y and the matrix L, the output is the row-vector

end