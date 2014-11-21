function G = findG(X,Y,K,L,rBF,params)
%
% findG Numerically evaluate Abel Transform basis functions used in the pBASEX
% inversion algorithm.
%
% G = findG(X,Y,K,L,params) calculates the value of the Abel transform of
% the numel(K)*numel(L) basis functions at each data pixel (x,y).
%
% Inputs:
%
% X and Y are 1-D arrays specifying the cartesian grid on which to sample
% the Abel transform of the basis functions.
%
% K is a 1-D array indexing the radial basis functions.
%
% L is a 1-D array of non-negative and even integers indexing the angular
% basis functions, which are the Legendre polynomials.
%
% rBF is a handle to a function that takes as inputs a number x,a number k,
% and a 1-D array params and outputs the value of the radial basis function
% specified by k and params at the point x.
%
% params is a 1-D array of parameters of the radial basis functions.
%
% Output:
%
% G is a 2-D array with the point g(a,b) being the value of the Abel transform
% of the basis function fkl(R,theta) = rBF(R,k,params)*legendre(l,m=0,th)
% at the point (x,y) with k = K(1+floor(b/numel(L))), l =
% L(1+mod(b,numel(L))), x = X(1+mod(a,numel(Y))), and y =
% Y(1+floor(a,numel(Y))).
%
% Example:
%
% % 256x256 pixels of data
% X = 0:255;
% Y = 0:255;
% % One basis function each 2 pixels
% K = 1:2:255;
% % Two-photon process can be represented angularly by the Legendre
% % polynomials with l = 0,2,4.
% L = 0:2:4;
% % Gaussian radial basis functions centered at k and of width sigma =
% % params(1).
% rBF = @(x,k,params) exp(-(x-k).^2/(2*params(1)^2))/k^2;
% params = sqrt(2);
% G = findG(X,Y,K,L,rBF,params);

zIP = @(r,k,params) sqrt((sqrt(2*10)*params(1)+k)^2-r^2); % Zero integrand point

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

u = 0:0.1:zIP(min(R),max(K),params); % Points at which to sample the integrand for the trapz integration
lenU = numel(u);

% Set up progress bar
progStep = ceil(lenX*lenY/500)+1;
progBar = ParforProgMon('Polar Integrals Progress:', lenX*lenY, progStep, 400, 70);

G = zeros(lenX*lenY,lenK*lenL); % Initialize output matrix

% Calculate integrals
parfor ind = 1:lenX*lenY % Loop over every data pixel
    
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
    
    if not(mod(ind,progStep))
        progBar.increment();
    end
    
end

progBar.delete();

end

function p=leg(l,x)

% leg Calculate the value of the l-th order Legendre polynomial (m=0) at
% point x. Implemented only for even values of l.
%
% p = leg(l,x)

switch l
    case 0
        p=ones(size(x));
        return
    case 1
        p=x;
        return
    case 2
        p=(3*x.*x -1)/2;
        return
    case 4
        x2=x.*x;
        p = ((35.*x2-30).*x2+3)/8;
        return
    case 6
        x2=x.*x;
        p = (((231.*x2-315).*x2+105).*x2-5)/16;
        return
    case 8
        x2=x.*x;
        p = ((((6435.*x2-12012).*x2+6930).*x2-1260).*x2+35)/128;
        return
    case 10
        x2=x.*x;
        p = (((((46189.*x2-109395).*x2+90090).*x2-30030).*x2+3465).*x2-63)/256;
        return
    case 12
        x2=x.*x;
        p = ((((((676039.*x2-1939938).*x2+2078505).*x2-1021020).*x2+225225).*x2-18018).*x2+231)/1024;
        return
    case 14
        x2=x.*x;
        p = (((((((5014575.*x2-16900975).*x2+22309287).*x2-14549535).*x2+4849845).*x2-765765).*x2+45045).*x2-429)/2048;
        return
    case 16
        x2=x.*x;
        p = ((((((((300540195.*x2-1163381400).*x2+1825305300).*x2-1487285800).*x2+669278610).*x2-162954792).*x2+19399380).*x2-875160).*x2+6435)/32768;
        return
    case 18
        x2=x.*x;
        p = (((((((((2268783825.*x2-9917826435).*x2+18032411700).*x2-17644617900).*x2+10039179150).*x2-3346393050).*x2+624660036).*x2-58198140).*x2+2078505).*x2-12155)/65536;
        return
    case 20
        x2=x.*x;
        p = ((((((((((34461632205.*x2-167890003050).*x2+347123925225).*x2-396713057400).*x2+273491577450).*x2-116454478140).*x2+30117537450).*x2-4461857400).*x2+334639305).*x2-9699690).*x2+46189)/262144;
        return
    case 22
        x2=x.*x;
        p = (((((((((((263012370465.*x2-1412926920405).*x2+3273855059475).*x2-4281195077775).*x2+3471239252250).*x2-1805044411170).*x2+601681470390).*x2-124772655150).*x2+15058768725).*x2-929553625).*x2+22309287).*x2-88179)/524288;
        return
    case 24
        x2=x.*x;
        p = ((((((((((((8061900920775.*x2-47342226683700).*x2+121511715154830).*x2-178970743251300).*x2+166966608033225).*x2-102748681866600).*x2+42117702927300).*x2-11345993441640).*x2+1933976154825).*x2-194090796900).*x2+10039179150).*x2-202811700).*x2+676039)/4194304;
        return
end
end