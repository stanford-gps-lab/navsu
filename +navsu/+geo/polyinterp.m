function [yy, yy_dot, chi2, p] = polyinterp(x, y, m_order, xx, flag, var, polyIn, nearestAdjust)

% SYNTAX:
%   [yy, yy_dot, chi2, p] = polyinterp(x, y, ord, xx, flag, var, polyIn, nearestAdjust)
%
% INPUT:
%   x   n x 1 <double> vector of x-data values
%   y   n x m <double> matrix of y-data values
%   ord 1 x 1 <int> order of the polynomial to be fit
%   xx  p x 1 <double> vector of evaluation points
% Optional inputs:
%   flag    n x 1 <logical> vector that is true for data points to be
%           ignored in this fit. Default = false
%   var     n x m <double> matrix of e.g. measurement variances. This 
%           scales the chi^2 statistic output. Default = 1
%   polyIn  1 x ord+1 <double> pre-fit polynomial of order ord to be used
%           instead of a fit computed on the provided data. Default = []
%   nearestAdjust   1 x 1 <logical> if true the fit error of the nearest 
%                   neighbor data point affects the result. Default = false

% OUTPUT:
%   yy      p x m <double> matrix of interpolated data points
%   yy_dot  p x m <double> matrix of gradient at fit points
%   chi2    chi-squared statistic of fit error
%   p       exponents of polynomial used
%
% DESCRIPTION:
%   Polynomial interpolation algorithm.



if nargin < 6
    var = ones(size(y));
elseif all(var == 0) || isempty(var)
    var = ones(size(y));
end

idx  = ~isnan(x) & ~any(isnan(y), 2);
if nargin > 4
    %remove flagged data and NaN
    idx  = idx & ~flag;
end

if nargin < 7
   polyIn = []; 
end

if nargin < 8
    nearestAdjust = true;
end

x = x(idx);
y = y(idx, :);
var = var(idx);

% make sure enough points remain to fit
if size(y, 1) < m_order + 1
%     warning('Warning bad fit in polyinterp\n');
    yy     = NaN(length(xx), size(y, 2));
    yy_dot = NaN(length(xx), size(y, 2));
    chi2   = NaN;    
    p = NaN(size(y, 2), m_order+1);
    return
end



xscale = range(x) / 2;
xbias  = mean(x);
xp = (x-xbias) / xscale;

yscale = range(y, 1) / 2;
ybias  = mean(y, 1);
if yscale ~= 0
    yp = (y-ybias) ./ yscale;
else
    yp = y*0;
end

n = size(xx, 1);
xx0 = xx;
xx = [xx; x];

if isempty(polyIn)
    p = navsu.geo.polyfit2(xp, yp, m_order);
else
    p = polyIn;
end
z = yscale .* (((xx - xbias)/xscale).^(m_order:-1:0) * p') + ybias; % vectorized version polyval
yy = z(1:n, :);

if nearestAdjust
    [~, indNearest] = min(abs(xx0 - x));
    
    yAdjust = z(n+indNearest, :) - y(indNearest, :);
    yy = yy - yAdjust;
end

if nargout > 1
    pdot = p(:, 1:end-1).*(m_order:-1:1);
    
    yy_dot = yscale .* (((xx(1:n) - xbias)/xscale).^(m_order-1:-1:0) * pdot') / xscale;
%     yy_dot = yscale*polyval(pdot, (xx(1:n) - xbias)/xscale)/xscale;
    
    if nargout > 2
        chi2 = dot(y - z((n+1):end, :), (y - z((n+1):end, :))./var, 2);
    end
end


