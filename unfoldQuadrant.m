function M_out = unfoldQuadrant(M,center)
%
% unfoldQuadrant Generate a full image from one quadrant.
%
% M_out = unfoldQuadrant(M, center) returns the full image with up-down and
% left-right symmetry given the lower-right quadrant.
%
% Inputs:
%
% M, a 2-D array representing the fourth quadrant of a symmetric image.
%
% center, a boolean that controls where the center of the image is. If
% true, the center of the first entry is the center of the image. If false,
% the upper-left corner of the first entry is the center of the image. If
% not specified, center defaults to True.
%
% Outputs:
%
% M_out, a 2-D array representing the full image.

if nargin<2
    center = 1;
end

M_out = [M(end:-1:1+center,end:-1:1+center,:),M(end:-1:1+center,:,:);M(:,end:-1:1+center,:),M]./4;

end