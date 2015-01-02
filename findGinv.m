function [Ginv,gData] = findGinv(gData, progBar)
%
% findGinv Numerical evaluation of polar basis functions in a
% cartesian basis.
%
% [Ginv, gData] = findGinv(gData, progBar) calculates the value of the a
% set of basis functions at specified image pixels.
%
% Inputs:
%
% gData, a structure with fields:
%
% *x and y, 1-D arrays specifying the cartesian grid on which to sample the
% Abel transform of the basis functions.
%
% *k, a 1-D array indexing the radial basis functions.
%
% *l, a 1-D array of non-negative and even integers indexing the angular
% basis functions, which are the Legendre polynomials.
%
% *rBF, a handle to a function that takes as inputs a number r, a number k,
% and an object params and outputs the value of the radial basis function
% specified by k and params at the radius r.
%
% *params, an object of any type used to specify parameters used in the
% radial basis functions.
%
% progBar, a boolean that controls whether or not to display a progress
% bar. If not specified, the default is False.
%
% Output:
%
% Ginv, a 2-D array with the point g(a,b) being the value of the basis
% function fkl(R,theta) = rBF(R,K,params)*legendre(L,m=0,th) at the point
% (X,Y) with K = k(1+floor(b/numel(l))), L = l(1+mod(b,numel(l))), X =
% x(1+mod(a,numel(y))), and Y = y(1+floor(a,numel(y))).
%
% gData, the same structure as in the input, unless some of the fields were
% missing in which case this output structure is updated with the defaults
% used for the calculation.
%
% Example:
%
% % 512x512 pixels of data
% gData.x = 0:511; gData.y = 0:511;
% % One basis function each 2 pixels
% gData.k = 0.5:2:511;
% % Two-photon process can be represented angularly by the Legendre
% % polynomials with l = 0,2,4
% gData.l = 0:2:4;
% % Gaussian radial basis functions centered at k and of width sqrt(2)
% % pixels
% gData.params = sqrt(2);
% gData.rBF = @(x,k,params) exp(-(x-k).^2/(2*params(1)^2))/k^2;
% gData.zIP = @(r,k,params) sqrt((sqrt(2*10)*params(1)+k).^2-r^2);
% Ginv = findGinv(gData);

% Set unspecified inputs to their default values
if nargin==1
    progBar = 0;
end
if ~isfield(gData,'x')
    gData.x = 0:255;
end
if ~isfield(gData,'y')
    gData.y = gData.x;
end
if ~isfield(gData,'k')
    gData.k = (gData.x(1:2:end)+gData.x(2:2:end))/2;
end
if ~isfield(gData,'l')
    gData.l = [0,2,4];
end
if ~isfield(gData,'rBF')
    gData.rBF = @(x,k,params) exp(-(x-k).^2/(2*params(1)^2))/k^2;
    gData.params = 2;
    gData.zIP = @(r,k,params) sqrt((sqrt(2*10)*params(1)+k).^2-r^2);
end
if ~isfield(gData,'params')
    gData.params = 1;
end

% Deal out gData values from structure to separate variables for code
% clarity and small (possibly negligible) communication decrease in parfor
X = gData.x;
Y = gData.y;
K = gData.k;
L = gData.l;
rBF = gData.rBF;
params = gData.params;

if nargin==1
    progBar = 0;
end

% Size of inputs
lenX = numel(X);
lenY = numel(Y);
lenK = numel(K);
lenL = numel(L);

% Create  (lenX*lenY)x1 arrays with the x, y, and R values for each data
% pixel
[xl,yl] = meshgrid(X,Y);
xl = xl(:);
yl = yl(:);
R = sqrt(xl.^2+yl.^2);
costh = yl./R;

% Set up progress bar
if progBar
    progStep = ceil(lenK*lenL/500)+1;
    progBar = ParforProgMon('Polar Integrals Progress:', lenK*lenL, progStep, 400, 70);
end

Ginv = zeros(lenK*lenL,lenX*lenY); % Initialize output matrix

% Calculate integrals
parfor ind = 1:lenK*lenL % Loop over every radial basis function
    
    k = K(ceil(ind/lenL));
    l = L(mod(ind-1,lenL)+1);
    
    Ginv(ind,:) = rBF(R,k,params).*leg(l,costh); % Calculate function value
    
    % Update progress bar
    if progBar&&not(mod(ind,progStep))
        progBar.increment();
    end
    
end

% Delete progress bar
if progBar
    progBar.delete();
end

Ginv(isnan(Ginv)) = 0;
Ginv = transpose(Ginv);

end