function [G,gData] = findG(gData, progBar)
%
% findG Numerical evaluation of Abel transformed polar basis functions in a
% cartesian basis.
%
% [G, gData] = findG(gData, progBar) calculates the value of the Abel
% transform of a set of basis functions at specified image pixels. The
% necessary integrals are done numerically using the trapezoid method.
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
% G, a 2-D array with the point g(a,b) being the value of the Abel
% transform of the basis function fkl(R,theta) =
% rBF(R,K,params)*legendre(L,m=0,th) at the point (X,Y) with K =
% k(1+floor(b/numel(l))), L = l(1+mod(b,numel(l))), X =
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
% G = findG(gData);

% Set unspecified inputs to their default values
if nargin==1
    progBar = 0;
end
if ~isfield(gData,'x')
    gData.x = 0:511;
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
    gData.params = sqrt(2);
end
if ~isfield(gData,'zIP')
    gData.zIP = @(r,k,params) 2*max(gData.x);
end
if ~isfield(gData,'trapzStep')
    gData.trapzStep = min([diff(gData.x),diff(gData.y)])/10;
    if ~gData.trapzStep
        gData.trapzStep = 0.1;
    end
end

% Deal out gData values from structure to separate variables for code
% clarity and small (possibly negligible) communication decrease in parfor
X = gData.x;
Y = gData.y;
K = gData.k;
L = gData.l;
rBF = gData.rBF;
params = gData.params;
zIP = gData.zIP;
trapzStep = gData.trapzStep;

% Size of inputs
lenX = numel(X);
lenY = numel(Y);
lenK = numel(K);
lenL = numel(L);

% Find radius of each pixel
[xl,yl] = meshgrid(X,Y);
xl = xl(:);
yl = yl(:);
R = sqrt(xl.^2+yl.^2);

% Points at which to sample the integrand for the trapz integration
u = 0:trapzStep:zIP(min(R),max(K),params);
lenU = numel(u);

% Set up progress bar
if progBar
    progStep = ceil(lenX*lenY/500)+1;
<<<<<<< HEAD
    pB = ParforProgMon('Polar Integrals Progress:', lenX*lenY, progStep, 400, 70);
else
    progStep = 0;
    pB = 0;
=======
    progBar = ParforProgMon('Polar Integrals Progress:', lenX*lenY, progStep, 400, 70);
else
    progStep = 0;
>>>>>>> 629229caabfed460d2a5fd6b8e5d81782bcc3051
end

G = zeros(lenX*lenY,lenK*lenL); % Initialize output matrix

% Calculate integrals
for ind = 1:lenX*lenY % Loop over every data pixel
    
    subG = zeros(1,lenK*lenL); % Initialize G subarray for for pixel (x,y)
    
    % Angular contribution to the integrand values (independent of k)
    cos_term = yl(ind)./sqrt(u.^2+R(ind)^2);
    cos_term(isnan(cos_term)) = 0; % cos(theta) = 0/0 at origin
    leg_terms = zeros(lenL,lenU);
    for lind = 1:numel(L)
        leg_terms(lind,:) = leg(L(lind),cos_term);
    end
    
    for q = 1:lenK % Loop over every radial basis function
        
        subU = u(u<zIP(R(ind),K(q),params)); % Shorten integration range
        
        if ~isempty(subU)
            
            rad_term = rBF(sqrt(subU.^2+R(ind)^2),K(q),params); % Radial contribution to the integrand
            
            for lind = 1:numel(L) % Loop over every angular basis function
                subG((q-1)*lenL+lind) = trapz(subU,rad_term.*leg_terms(lind,1:numel(subU)),2); % Calculate the integral
            end
            
        end
        
    end
    
    G(ind,:) = subG;
    
    % Update progress bar
    if progBar&&not(mod(ind,progStep))
        pB.increment();
    end
    
end

% Delete progress bar
if progBar
    pB.delete();
end

end