function out = pbasex(images, gData, makeImages, weights, regularization, alpha)
%
% pbasex Apply a modification of the PBASEX algorithm (see Garcia et al.,
% Rev. Sci. Instrum. 75, 4989 (2004)) to Abel invert an image.
%
% out = pbasex(im,gData,makeImages,r) returns the inversion basis
% function coefficients and other inversion information for the image im
% using data in gData to perform the inversion.
%
% Inputs:
%
% im is a 2-D matrix representing one quadrant of an image of Abel
% transformed data. A 3-D matrix can be passed as well, where each entry in
% the third dimension is an individual image.
%
% gData is either a filename pointing to a .mat file with the relevant
% inversion data or a structure object with the relevant inversion data as
% fields. The inversion data fields are:
%
%   (1) x and y, 1-D arrays where the i-th element is the distance from the
% image center to the i-th pixel in each dimension. Dimensions must match
% such that size(im,1) = size(x) and size(im,2) = size(y).
%
%   (2) k, a 1-D array indexing the radial part of the basis functions.
%
%   (3) l, a 1-D array indexing the angular part of the basis functions.
%   The
% values of l must be non-negative, even integers.
%
%   (4) rBF, a function handle where rBF(r,k,params) outputs the value of
% the k-th radial basis function at radius r.
%
%   (5) params, an object holding the parameters to pass to rBF.
%
%   (6) Up and V, 2-D matrices, and Sinv, a 1-D array, holding the data of
% the Singular Value Decomposition of the 2-D matrix G such that G =
% Up'*diag(S)*V'. G holds the values of the Abel transform of all the basis
% functions at each image pixel. G(a,b) is the value of the Abel transform
% of the basis function fkl(R,theta) =
% rBF(R,kval,params)*legendre(lval,m=0,th) at the point (xval,yval) with
% kval = k(1+floor(b/numel(l))), lval = l(1+mod(b,numel(l))) , xval =
% x(1+mod(a,numel(y))), and yval = y(1+floor(a,numel(y))). See findG.m for
% more details.
%
%   (7) Ginv, a 2-D matrix holding the values of all the basis functions at
% each image pixel. Ginv(a,b) is the value of the basis function
% fkl(R,theta) = rBF(R,kval,params)*legendre(lval,m=0,th) at the point
% (xval,yval) with kval = k(1+floor(b/numel(l))), lval =
% l(1+mod(b,numel(l))) , xval = x(1+mod(a,numel(y))), and yval =
% y(1+floor(a,numel(y))). Ginv is only needed if makeImages is True. See
% findGinv.m for more details.
%
% makeImages, a boolean that controls whether or not images for the fitted
% data and the fitted inverted data are created. If not specified, the
% default is to not make these images.
%
% r, a 1-D array of values at which to sample the fitted radial energy
% spectrum and beta values. If not specified, the default is to use
% gData.x.
%
% Outputs:
%
% out_struct is a structure holding the inversion information, with fields:
%
%   (1) r, a 1-D array holding the value, in pixels, of the radii at which
% the calculated radial energy spectrum was sampled.
%
%   (2) Ir, a 1-D array holding the value of the calculated radial energy
% spectrum at the radii specified in r. If multiple images were inverted,
% Ir becomes a 2-D matrix where each column is the calculated radial energy
% spectrum for an image.
%
%   (3) k, a 1-D array indexing the radial part of the basis functions.
%
%   (4) c, a 2-D matrix holding the inversion data. c(a,b) corresponds to
% the weight given to the basis function indexed in the radial coordinate
% by k(a) and in the polar coordinate by l(b). If multiple images were
% inverted, c becomes a 3-D matrix where the third dimension indexes the
% images.
%
%   (5) betas, a 2-D matrix holding the values of the beta parameters of
% the fit. betas(a,b) corresponds to the value of the beta parameter of
% order l(b+1) at the radial coordinate x(a). If multiple images were
% inverted, betas becomes a 3-D matrix where the third dimension indexes
% the images.
%
%   (6) recon and inv, 2-D matrices holding images of the reconstructed and
% inverted data. If multiple images were inverted, they become 3-D matrices
% where the third dimension indexes the images. These are only created if
% makeImages is True.

% By default, do not use regularization
if nargin<3
    makeImages = 0;
end

if nargin<4
    weights = 0;
end

% By default, do not make reconstruction images
if nargin<5
    regularization = 0;
end

% Set default alpha
if nargin<6
    alpha = 4e-5;
end

% Load gData if file specified
if ischar(gData)
    gData = loadG(gData);
end

% Problem Dimensionality
nx = numel(gData.x);
nk = gData.nk;
nl = gData.nl;
nims = size(images,3);

% Invert the data
images = reshape(images,nx^2,nims);
if any(weights)
    c = gData.V*(diag(gData.S./(gData.S.^2+regularization))*((gData.Up*bsxfun(@times,U,w(:)))\(gData.Up*bsxfun(@times,images,w(:)))));
else
    c = gData.V*(diag(gData.S./(gData.S.^2+regularization))*(gData.Up*images));
end

% Calculate the radial intensity and beta values
E = alpha*gData.x.^2;
IEB = 1/(2*alpha)*diag(gData.x)*(gData.frk*reshape(permute(reshape(c,nk,nl,nims),[1,3,2]),nk,nl*nims));
IE = IEB(:,1:nims);
betas = bsxfun(@times,permute(reshape(IEB(:,nims+1:end),[nx,nims,nl-1]),[1,3,2]),1./permute(IE,[1,3,2]));

% Generate Abel transformed and phi=0 sliced images from fit
if makeImages
    im_recon = unfoldQuadrant(reshape(gData.Up'*(diag(1./gData.Sinv)*(gData.V'*c)),nx,nx,nims));
    im_inv = unfoldQuadrant(reshape(gData.Ginv*c,nx,nx,nims));
end

% Make an output structure with commonly used data
out = struct('E',E,'IE',IE,'betas',betas,'c',c);
if makeImages
    out.recon = im_recon;
    out.inv = im_inv;
end

end