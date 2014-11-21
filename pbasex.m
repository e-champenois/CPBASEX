function out_struct = pbasex(im,gData)
%
% pbasex Apply the PBASEX algorithm (see ref. [1]) to Abel invert an image.
%
% out_struct = pbasex(im,gData) returns the inversion basis function coefficients
% and other inversion information for the image im using data in gData to
% perform the inversion.
%
% Inputs:
%
% im is a 2-D array representing one quadrant of an image of Abel
% transformed data. A 3-D array can be passed as well, where each entry in
% the third dimension is an individual image.
%
% gData is either a filename pointing to a .mat file with the relevant
% inversion data or a structure object with the relevant inversion data as
% fields. The inversion data fields are:
%
% (1) x and y, 1-D arrays where the i-th element is the distance from the
% image center to the i-th pixel in each dimension. Dimensions must match
% such that size(im,1) = size(x) and size(im,2) = size(y).
%
% (2) k, 1-D array indexing the radial part of the basis functions.
%
% (3) l, 1-D array indexing the angular part of the basis functions. The
% values of l must be non-negative, even integers.
%
% (4) rBF, a function handle where rBF(r,k,params) outputs the value of the
% k-th radial basis function at radius r.
%
% (5) params, parameters to pass to rBF.
%
% (6) G, a 2-D array holding the Abel transform of every basis function.
% G(a,b) is the value of the Abel transform of the basis function
% fkl(R,theta) = rBF(R,kval,params)*legendre(lval,m=0,th) at the point
% (xval,yval) with kval = k(1+floor(b/numel(l))), lval =
% l(1+mod(b,numel(l))) , xval = x(1+mod(a,numel(y))), and yval =
% y(1+floor(a,numel(y))). See findG.m for more details.
%
% (7) Ginv, a 2-D array holding all the basis functions. Ginv(a,b) is the
% value of the basis function fkl(R,theta) =
% rBF(R,kval,params)*legendre(lval,m=0,th) at the point (xval,yval) with
% kval = k(1+floor(b/numel(l))), lval = l(1+mod(b,numel(l))) , xval =
% x(1+mod(a,numel(y))), and yval = y(1+floor(a,numel(y))). See findG.m for
% more details.
%
% Outputs:
%
% out_struct is a structure holding the inversion information, with fields:
%
% (1) r, a 1-D array of 

% Load gData if file specified
if ischar(gData)
    gData = load(gData);
end

% Problem dimensions
lenX = numel(gData.x);
lenY = numel(gData.y);
lenK = numel(gData.k);
lenL = numel(gData.l);
numIms = size(im,3);

% Sample radial part of basis functions: f(r,k)
[K,Rm]=meshgrid(gData.k,gData.x);
fk = gData.rBF(Rm,K,gData.params);

% Preallocate memory for C, I_R, beta, im_recon, im_inv
C = zeros(lenK,lenL,numIms);
I_R = zeros(lenX,numIms);
beta = zeros(lenX,lenL-1,numIms);
im_recon = zeros(2*lenX-1,2*lenX-1,numIms);
im_inv = im_recon;

% Invert the data one image at a time
for i = 1:numIms % Loop over every image
    
    % Find basis function coefficients
    im1d = reshape(im(:,:,i),lenX*lenY,1);
    c = gData.V*(gData.Up*im1d./gData.S);
    C(:,:,i)=reshape(c,lenL,lenK)';
    
    % Calculate radial energy spectrum
    I_R(:,i) = gData.x.^2.*(fk*C(:,1,i))';
    
    % Calculate beta parameters
    for lval = 1:lenL-1
        beta(:,lval,i) = gData.x'.^2.*(fk*C(:,lval+1,i))./I_R(:,i);
    end
    
    % Generate Abel transformed and phi=0 sliced images from fit
    im_recon(:,:,i) = unfoldQuadrant(reshape(gData.G*c,lenX,lenY));
    im_inv(:,:,i) = unfoldQuadrant(reshape(gData.Ginv*c,lenX,lenY));
    
end

% Make an output structure with commonly used data
out_struct = struct('r',gData.x,'Ir',I_R,'k',gData.k,'c',C,'beta',beta,'recon',im_recon,'inv',im_inv);

end
